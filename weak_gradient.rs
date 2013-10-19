use common::*;
use vector_monomial::VectorMonomial;
use monomial::{Monomial, DegLim, MaxMonDeg, MaxMonFactorDeg};
use polynomial::{PolyBorrowing};
use mesh::{Mesh, OShape, SideFace};
use dense_matrix::DenseMatrix;
use lapack;
use lapack::lapack_int;

use std::vec;

/*
 * For a weak function v on a finite element T, the weak gradient of degree r
 * of v on T is defined to be the polynomial vector function wgrad(v) in
 * [P_r(T)]^d, such that:
 *
 * [WGRAD_DEF]
 *   (wgrad(v), q)_T = -(v_0, div q)_T + <v_b, q.n>_bnd(T), for all q in [P_r(T)]^d
 *
 * By linearity of both sides in the right sides of the inner products, the above
 * holds iff the same equation holds with the q's further restricted to be
 * vector monomials on T (which are a monomial in one range component and 0 for
 * all other components), which functions form a basis {b_i}_i for [P_r(T)]^d.
 *
 * Letting q = b_i, and writing wgrad(v) = sum_j eta_j b_j, we obtain a linear
 * system from WGRAD_DEF which we can solve for the unkowns eta_j, with
 * (b_i, b_j)_T being the matrix elements of the system, and the right hand side
 * of WGRAD_DEF defining the right hand side column vector for the system.
 *
 * We will only actually need weak gradients of weak functions which are monomials
 * on a single face (interior or side) of a finite element and 0 elsewhere. These
 * are given by the functions wgrad_int_mon and wgrad_side_mon.
 */

pub struct WeakGrad {
  comp_mon_coefs: ~[~[R]]
}

pub struct WeakGradSolver<Mon> {

  wgrad_comp_mons_deg_lim: DegLim,
  
  wgrad_comp_mons: ~[Mon],

  basis_vmons: ~[VectorMonomial<Mon>], // ordered ascending by (monomial component dimension, monomial exponents tuple).

  ips_basis_vmons_by_oshape: ~[DenseMatrix], // Basis inner products are in upper triangular parts of the matrices.

  // Work matrices for lapack to avoid allocations and because lapack enjoys writing over its inputs.
  lapack_ips_basis_vmons: ~DenseMatrix,
  lapack_pivots: ~[lapack_int],
  lapack_pivots_buf: *mut lapack_int,
  lapack_rhs: ~DenseMatrix,
}

impl <Mon:Monomial> WeakGradSolver<Mon> {

  pub fn new<MESHT:Mesh<Mon>>(comp_mons_deg_lim: DegLim, mesh: &MESHT) -> WeakGradSolver<Mon> {
    lapack::init(); // TODO: Move this to main when available.
    let comp_mons: ~[Mon] = Monomial::mons_with_deg_lim_asc(comp_mons_deg_lim);
    let vmons = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(comp_mons);
    let num_vmons = vmons.len();
    
    let ips = vec::from_fn(mesh.num_oriented_element_shapes(), |os| {
      DenseMatrix::upper_triangle_from_fn(num_vmons, |i,j| {
        if vmons[i].mon_dim != vmons[j].mon_dim { 0 as R }
        else {
          mesh.intg_facerel_mon_on_oshape_int(vmons[i].mon * vmons[j].mon, OShape(os))
        }
      })
    });

    let mut lapack_pivots = vec::from_elem(num_vmons, 0 as lapack_int);
    let lapack_pivots_buf = vec::raw::to_mut_ptr(lapack_pivots);

    WeakGradSolver {
      wgrad_comp_mons_deg_lim: comp_mons_deg_lim,
      wgrad_comp_mons: comp_mons,
      basis_vmons: vmons,
      ips_basis_vmons_by_oshape: ips,
      lapack_ips_basis_vmons: ~DenseMatrix::from_elem(num_vmons, num_vmons, 0 as R),
      lapack_pivots: lapack_pivots,
      lapack_pivots_buf: lapack_pivots_buf,
      lapack_rhs: ~DenseMatrix::from_elem(num_vmons, 200, 0 as R) // Initially allocate for up to 200 shape funs per oshape -
                                                                  // will reallocate if necessary.
    }
  }

 /*
  * These two functions compute one component of the right hand side of the equation (WGRAD_DEF),
  *   WGRAD_DEF_RHS:    -(v_0, div q)_T + <v_b, q.n>_bnd(T),
  * on a reference finite element of the mesh for an interior or side supported monomial v.
  */
  fn wgrad_def_rhs_for_int_mon<MESHT:Mesh<Mon>>(&self, v: Mon, oshape: OShape, q: &VectorMonomial<Mon>, mesh: &MESHT) -> R {
    // Interior supported v: only the -(v_0, div q)_T term can be non-zero in the rhs of (WGRAD_DEF).
    let (div_q_coef, div_q_mon) = q.divergence_coef_and_mon();
    -div_q_coef * mesh.intg_facerel_mon_on_oshape_int(v * div_q_mon, oshape)
  }
  
  fn wgrad_def_rhs_for_side_mon<MESHT:Mesh<Mon>>(&self, v: Mon, oshape: OShape, side_face: SideFace,
                                                        q: &VectorMonomial<Mon>, mesh: &MESHT) -> R {
    // Side supported v: only the <v_b, q.n>_bnd(T) term can be non-zero in the rhs of (WGRAD_DEF).
    mesh.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(v, q, oshape, side_face)
  }

  #[fixed_stack_segment] 
  #[inline(never)]
  pub fn wgrads_on_oshape<MESHT:Mesh<Mon>>(&mut self, int_mons: &[Mon], side_mons_by_side: &[&[Mon]],
                                                      oshape: OShape, mesh: &MESHT) -> (~[WeakGrad], ~[~[WeakGrad]]) {
    let num_vmons = self.basis_vmons.len();

    let sols_col_maj = 
      unsafe {
        // a - The system matrix of vmon inner products [in] linearized in column-major order.
        let a = {
          self.ips_basis_vmons_by_oshape[*oshape].copy_upper_triangle_into(self.lapack_ips_basis_vmons); 
          self.lapack_ips_basis_vmons.mut_col_maj_data_ptr()
        };
        
        // b - The WGRAD_DEF system right hand side column vectors matrix [in], solution columns [out], both
        //     linearized in column-major order.
        let (b, num_rhs_cols) = {
          let rhss = self.wgrad_def_rhss(int_mons, side_mons_by_side, oshape, mesh);
          (rhss.mut_col_maj_data_ptr(), rhss.num_cols)
        };

        lapack::solve_symmetric_as_col_maj_with_ut_sys(a, num_vmons as lapack_int,
                                                       b, num_rhs_cols as lapack_int,
                                                       self.lapack_pivots_buf);
        
        vec::from_buf(b as *R, num_vmons * num_rhs_cols)
      };

    // Interior weak gradients are represented as a vector, indexed by interior shape function (monomial) number,
    // of polynomial vectors representing the weak gradient of the shape function.
    let int_wgrads = self.int_wgrads_from_combined_sol_coefs(sols_col_maj.slice(0, num_vmons * int_mons.len()));

    // Side weak gradients are represented as a vector, indexed by side number, of vectors, indexed by side shape function
    // (monomial) number, of polynomial vectors representing the weak gradient of the shape function.
    let num_side_mons_by_side = side_mons_by_side.map(|mons| mons.len());
    let side_wgrads_by_side = self.side_wgrads_from_combined_sol_coefs(sols_col_maj.slice_from(num_vmons*int_mons.len()),
                                                                       num_side_mons_by_side);

    (int_wgrads, side_wgrads_by_side)
  }
 
  // Pack WGRAD_DEF right hand side computations as column vectors in a combined rhs matrix suitable for an LAPACK solver.
  // Returns a matrix consisting of right hand side column vectors of the WGRAD_DEF system evaluated over the sequence of
  // basis vector monomials which are represented by the rows.  The first columns represent the interior supported
  // monomials on the oshape, ordered by increasing monomial. These are followed by columns representing the monomials
  // supported on the sides in order of increasing side number, and by increasing monomial within a side section.
  fn wgrad_def_rhss<'a,MESHT:Mesh<Mon>>(&'a mut self, int_mons: &[Mon], side_mons_by_side: &[&[Mon]],
                                                      oshape: OShape, mesh: &MESHT) -> &'a mut DenseMatrix {
    let num_vmons = self.basis_vmons.len();
    let num_rhs_cols = {
      let total_side_mons_all_sides = side_mons_by_side.iter().fold(0u, |sum, side_mons| sum + side_mons.len());
      int_mons.len() + total_side_mons_all_sides
    };

    if self.lapack_rhs.capacity_cols < num_rhs_cols {
      self.lapack_rhs = ~DenseMatrix::of_size_with_cols_capacity(num_vmons, num_rhs_cols, 2*num_rhs_cols);
    } else {
      self.lapack_rhs.set_num_cols(num_rhs_cols);
    }

    for c in range(0, int_mons.len()) {
      let mon = int_mons[c].clone();
      for r in range(0, num_vmons) {
        let rhs_val = self.wgrad_def_rhs_for_int_mon(mon.clone(), oshape, &self.basis_vmons[r], mesh);
        self.lapack_rhs.set(r, c, rhs_val);
      }
    }
    let mut c = int_mons.len();
    for side_num in range(0, side_mons_by_side.len()) {
      let side_face = SideFace(side_num);
      let side_mons = side_mons_by_side[side_num];
      // Write a column of rhs data for each monomial on this side.
      for mon_num in range(0, side_mons.len()) {
        let mon = side_mons[mon_num].clone();
        for r in range(0, num_vmons) {
          let rhs_val = self.wgrad_def_rhs_for_side_mon(mon.clone(), oshape, side_face, &self.basis_vmons[r], mesh);
          self.lapack_rhs.set(r, c, rhs_val);
        }
        c += 1;
      }
    }
    
    &mut *self.lapack_rhs
  }
 
  // Unpack slice of solution coefficients for interior monomials from LAPACK solver as weak gradients.
  fn int_wgrads_from_combined_sol_coefs(&self, int_wgrad_coefs: &[R]) -> ~[WeakGrad] {
    int_wgrad_coefs
      .chunk_iter(self.basis_vmons.len()) // Chunk into sections corresponding to wgrad coefs of individual int mons.
      .map(|wgrad_vmon_coefs| 
           WeakGrad {
             comp_mon_coefs: wgrad_vmon_coefs.chunk_iter(self.wgrad_comp_mons.len()) // comp dim sections
                                             .map(|comp_coefs| comp_coefs.to_owned())
                                             .collect()
           })
      .collect()
  }
 
  // Unpack slice of solution coefficients for side monomials from LAPACK solver as weak gradients.
  fn side_wgrads_from_combined_sol_coefs(&self, sides_wgrad_coefs: &[R], num_side_mons_by_side: &[uint]) -> ~[~[WeakGrad]] {
    let num_vmons = self.basis_vmons.len();
    let side_start_ixs = cumulative_sums_prev_elems(num_side_mons_by_side.map(|&num_mons| num_vmons * num_mons));
    range(0, num_side_mons_by_side.len()).map(|side_num| {
      // Produce a vector of wgrads for this side.
      let side_start_ix = side_start_ixs[side_num];
      let side_wgrad_coefs = sides_wgrad_coefs.slice(side_start_ix, side_start_ix + num_vmons * num_side_mons_by_side[side_num]);
      side_wgrad_coefs           // Solution coefficients for all monomials on this side.
        .chunk_iter(num_vmons)   // Chunk into sections corresponding to wgrads of individual side monomials.
        .map(|wgrad_vmon_coefs|
             WeakGrad {
               comp_mon_coefs: wgrad_vmon_coefs.chunk_iter(self.wgrad_comp_mons.len()) // comp dim sections
                                               .map(|comp_coefs| comp_coefs.to_owned())
                                               .collect() // Collect component coefficients for a single weak grad.
             })
        .collect() // Collect vector of weak gradients by monomial number.
    }).collect()   // Collect weak gradient collections by side number.
  }


  pub fn weak_grad_ops(&self) -> ~WeakGradOps<Mon> {
    use std::hashmap::HashMap;

    let comp_mons = self.wgrad_comp_mons.clone(); // TWEAK: Could try borrowing this instead (add lifetime param to the type). 

    let prod_mons = { 
      let prod_mons_deg_lim = match self.wgrad_comp_mons_deg_lim {
        MaxMonDeg(l) => MaxMonDeg(2*l),
        MaxMonFactorDeg(l) => MaxMonFactorDeg(2*l) 
      };
      Monomial::mons_with_deg_lim_asc(prod_mons_deg_lim)
    };
      
    let (num_comp_mons, num_prod_mons) = (comp_mons.len(), prod_mons.len());

    let comp_monns_by_prod_monn = {
      let mut comp_monns_by_prod_monn: HashMap<Mon,~[(uint,uint)]> = HashMap::with_capacity(num_prod_mons);
      for fac1_monn in range(0, num_comp_mons) {
        for fac2_monn in range(fac1_monn, num_comp_mons) {
          let prod_mon = comp_mons[fac1_monn] * comp_mons[fac2_monn];
          let approx_twice_avg_num_fac_pairs = sq(num_comp_mons)/num_prod_mons;
          comp_monns_by_prod_monn.find_or_insert_with(prod_mon, |_| vec::with_capacity(approx_twice_avg_num_fac_pairs))
                                 .push((fac1_monn, fac2_monn));
        }
      }
      vec::from_fn(num_prod_mons, |prod_monn| { comp_monns_by_prod_monn.get(&prod_mons[prod_monn]).clone() })
    };

    ~WeakGradOps {
      dotprod_coefs_buf: vec::from_elem(num_prod_mons, 0 as R),
      dotprod_mons: prod_mons,
      wgrad_comp_mons: comp_mons,
      wgrad_comp_monns_by_dotprod_monn: comp_monns_by_prod_monn
    }
  }

} // WeakGradSolver impl

/// This class holds context information and work buffers used to perform efficient operations on weak gradients. 
pub struct WeakGradOps<Mon> {
  wgrad_comp_mons: ~[Mon],
  // dot product support
  dotprod_coefs_buf: ~[R],
  dotprod_mons: ~[Mon],
  wgrad_comp_monns_by_dotprod_monn: ~[~[(uint,uint)]]
}

impl<Mon:Monomial> WeakGradOps<Mon> {

  pub fn dot<'a>(&'a mut self, wgrad1: &WeakGrad, wgrad2: &WeakGrad) -> PolyBorrowing<'a,Mon> {
    let space_dims = wgrad1.comp_mon_coefs.len();
    assert!(space_dims == wgrad2.comp_mon_coefs.len());

    // Set the coefficients for the product monomials in the dot product coefficients buffer.
    for prod_monn in range(0, self.dotprod_mons.len()) {
      // Get the mon nums of the the component mons which produce this product monomial - only non-descending pairs included.
      let ref fac_monn_pairs = self.wgrad_comp_monns_by_dotprod_monn[prod_monn];
      // Sum the products of all coefficient pairs where the corresponding monomials multiply to give the product monomial.
      let sum_coef_prods = fac_monn_pairs.iter().fold(0 as R, |sum_coef_prods, &(fac1_monn, fac2_monn)|
        sum_coef_prods
          + range(0, space_dims).fold(0 as R, |dim_contrs, d| { // sum contributions for this factor pair over space dimensions
              dim_contrs
                + wgrad1.comp_mon_coefs[d][fac1_monn] * wgrad2.comp_mon_coefs[d][fac2_monn]
                + if fac1_monn != fac2_monn { wgrad1.comp_mon_coefs[d][fac2_monn] * wgrad2.comp_mon_coefs[d][fac1_monn] }
                  else { 0 as R }
            })
      );
      
      self.dotprod_coefs_buf[prod_monn] = sum_coef_prods;
    }
    
    PolyBorrowing::new(self.dotprod_coefs_buf, self.dotprod_mons)
  }
}


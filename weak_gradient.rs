use std::vec;
use std::vec::raw::to_ptr;

use common::*;
use vector_monomial::VectorMonomial;
use monomial::{Monomial, DegLim};
use mesh::{Mesh, OShape, SideFace};
use dense_matrix::DenseMatrix;
use std::libc::{c_double};

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

impl WeakGrad {

  pub fn lcomb(terms: &[(R,&WeakGrad)]) -> WeakGrad {
    if terms.len() == 0 { fail!("lcomb_wgrads: At least one weak gradient is required.") }
    let (space_dims, num_comp_mons) = match terms[0] { (_, wgrad) => (wgrad.comp_mon_coefs.len(), wgrad.comp_mon_coefs[0].len()) };
    WeakGrad {
      comp_mon_coefs: 
        vec::from_fn(space_dims, |d| 
          vec::from_fn(num_comp_mons, |mon_num|
            terms.iter().fold(0 as R, |sum, &(c, wgrad)| sum + c * wgrad.comp_mon_coefs[d][mon_num])))
    }
  }
  
}


pub struct OShapeWeakGrads {
  int_mon_wgrads: ~[WeakGrad],    // interior monomial weak gradients, indexed by interior monomial number
  side_mon_wgrads: ~[~[WeakGrad]] // sifn add(&self, rhs: &RHS) -> Result;de monomial weak gradients, indexed by side number then monomial number on the side
}


pub struct WeakGradSolver<M> {

  basis_vmons: ~[VectorMonomial<M>], // ordered ascending by (monomial component dimension, monomial exponents tuple).

  wgrad_comp_mons: ~[M],

  ips_basis_vmons_utmat_by_oshape: ~[DenseMatrix] // basis inner products are in upper triangular part of the matrix

}

impl <M:Monomial> WeakGradSolver<M> {

  pub fn new<MESHT:Mesh<M>>(deg_lim: DegLim, mesh: &MESHT) -> WeakGradSolver<M> {
    let comp_mons: ~[M] = Monomial::mons_with_deg_lim_asc(deg_lim);
    let vmons = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(comp_mons);
    let ips = vec::from_fn(mesh.num_oriented_element_shapes(), |os| {
      DenseMatrix::upper_triangular_from_fn(vmons.len(), vmons.len(), |i,j| {
        if vmons[i].mon_dim != vmons[j].mon_dim { 0 as R }
        else {
          mesh.intg_facerel_mon_on_oshape_int(vmons[i].mon * vmons[j].mon, OShape(os as u32))
        }
      })
    });
    WeakGradSolver { 
      basis_vmons: vmons,
      wgrad_comp_mons: comp_mons,
      ips_basis_vmons_utmat_by_oshape: ips
    }
  }

 /*
  * These two functions compute one component of the right hand side of the equation (WGRAD_DEF),
  *   WGRAD_DEF_RHS:    -(v_0, div q)_T + <v_b, q.n>_bnd(T),
  * on a reference finite element of the mesh for an interior or side supported monomial v.
  */
  fn wgrad_def_rhs_for_int_mon<M:Monomial,MESHT:Mesh<M>>(&self, v: M, oshape: OShape, q: &VectorMonomial<M>, mesh: &MESHT) -> R {
    // Interior supported v: only the -(v_0, div q)_T term can be non-zero in the rhs of (WGRAD_DEF).
    let (div_q_coef, div_q_mon) = q.divergence_coef_and_mon();
    -div_q_coef * mesh.intg_facerel_mon_on_oshape_int(v * div_q_mon, oshape)
  }
  fn wgrad_def_rhs_for_side_mon<M:Monomial,MESHT:Mesh<M>>(&self, v: M, oshape: OShape, side_face: SideFace, q: &VectorMonomial<M>, mesh: &MESHT) -> R {
    // Side supported v: only the <v_b, q.n>_bnd(T) term can be non-zero in the rhs of (WGRAD_DEF).
    mesh.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(v, q, oshape, side_face)
  }
  
  #[fixed_stack_segment] 
  #[inline(never)]
  pub fn wgrads_on_oshape<MESHT:Mesh<M>>(&self, int_mons: &[M], side_mons_by_side: &[&[M]], oshape: OShape, mesh: &MESHT) -> OShapeWeakGrads {

    // a - The system matrix of vmon inner products [in], factorization information [out], both linearized in column-major order.
    let lapack_a = self.ips_basis_vmons_utmat_by_oshape[*oshape].col_maj_data_copy();
    
    // b - The WGRAD_DEF system right hand side column vectors matrix [in], solution columns [out], both linearized in column-major order.
    let (lapack_b, nrhs) = self.wgrad_def_rhss_col_maj(int_mons, side_mons_by_side, oshape, mesh);
   
    let num_vmons = self.basis_vmons.len();
    let n = num_vmons as lapack_int;
    let lapack_ipiv = vec::from_elem(num_vmons, 0 as lapack_int); 
   
    unsafe {
      solve_symmetric_as_col_maj_with_ut_sys(to_ptr(lapack_a), n, to_ptr(lapack_b), nrhs as lapack_int, to_ptr(lapack_ipiv));
    }

    // Interior weak gradients are represented as a vector, indexed by interior shape function (monomial) number,
    // of polynomial vectors representing the weak gradient of the shape function.
    let int_wgrads = self.int_wgrads_from_combined_sol_coefs(lapack_b.slice(0, num_vmons * int_mons.len()));

    // Side weak gradients are represented as a vector, indexed by side number, of vectors, indexed by side shape function
    // (monomial) number, of polynomial vectors representing the weak gradient of the shape function.
    let num_side_mons_by_side = side_mons_by_side.map(|mons| mons.len());
    let side_wgrads_by_side = self.side_wgrads_from_combined_sol_coefs(lapack_b.slice_from(num_vmons * int_mons.len()), num_side_mons_by_side);

    OShapeWeakGrads {
      int_mon_wgrads: int_wgrads,
      side_mon_wgrads: side_wgrads_by_side
    }
  }
 
  // Pack WGRAD_DEF right hand side computations as column vectors in a combined rhs matrix suitable for an LAPACK solver.
  // Makes a linearized (column-major) matrix consisting of right hand side column vectors of the WGRAD_DEF system evaluated
  // over the sequence of basis vector monomials which are represented by the rows.  The first columns represent the
  // interior supported monomials on the oshape, ordered by increasing monomial. These are followed by columns representing
  // the monomials supported on the sides in order of increasing side number, and by increasing monomial within a side section.
  fn wgrad_def_rhss_col_maj<MESHT:Mesh<M>>(&self, int_mons: &[M], side_mons_by_side: &[&[M]], oshape: OShape, mesh: &MESHT) -> (~[R],uint) {
    let total_side_mons_all_sides = side_mons_by_side.iter().fold(0u, |sum, side_mons| sum + side_mons.len());
    let num_vmons = self.basis_vmons.len();
    let mut b = DenseMatrix::from_elem(num_vmons, int_mons.len() + total_side_mons_all_sides, 0 as R);
    for c in range(0, int_mons.len()) {
      let mon = int_mons[c].clone();
      for r in range(0, num_vmons) {
        b.set(r, c, self.wgrad_def_rhs_for_int_mon(mon.clone(), oshape, &self.basis_vmons[r], mesh));
      }
    }
    let mut c = int_mons.len();
    for side_num in range(0, side_mons_by_side.len()) {
      let side_face = SideFace(side_num as u8);
      let side_mons = side_mons_by_side[side_num];
      // Write a column of rhs data for each monomial on this side.
      for mon_num in range(0, side_mons.len()) {
        let mon = side_mons[mon_num].clone();
        for r in range(0, num_vmons) {
          b.set(r, c, self.wgrad_def_rhs_for_side_mon(mon.clone(), oshape, side_face, &self.basis_vmons[r], mesh));
        }
        c += 1;
      }
    }
    let num_rhss = b.num_cols;
    (b.data, num_rhss)
  }
 
  // Unpack slice of solution coefficients for interior monomials from LAPACK solver as weak gradients.
  fn int_wgrads_from_combined_sol_coefs(&self, int_wgrad_coefs: &[R]) -> ~[WeakGrad] {
    int_wgrad_coefs
      .chunk_iter(self.basis_vmons.len()) // Chunk into sections corresponding to wgrad coefs of individual interior monomials.
      .map(|wgrad_vmon_coefs| 
           WeakGrad {
             comp_mon_coefs: wgrad_vmon_coefs.chunk_iter(self.wgrad_comp_mons.len()) // Divide coefs into sections by component dimension.
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
      side_wgrad_coefs                                               // Solution coefficients for all monomials on this side.
        .chunk_iter(num_vmons)                                       // Chunk into sections corresponding to wgrads of individual side monomials.
        .map(|wgrad_vmon_coefs|
             WeakGrad {
               comp_mon_coefs: wgrad_vmon_coefs.chunk_iter(self.wgrad_comp_mons.len()) // Divide coefs into sections of common component dimension.
                                               .map(|comp_coefs| comp_coefs.to_owned())
                                               .collect() // Collect vector of component coefficients for a single weak gradient.
             })
        .collect() // Collect vector of weak gradients by monomial number.
    }).collect()   // Collect weak gradient collections by side number.
  }

} // WeakGradSolver impl

#[link_args = "-Llib/mkl lib/lapack_wrappers.o -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_core -lmkl_intel_thread -lmkl_core -liomp5"]
extern {

  fn solve_symmetric_as_col_maj_with_ut_sys(a: *c_double,
                                            n: lapack_int,
                                            b: *c_double,
                                            nrhs: lapack_int,
                                            ipiv: *lapack_int) -> lapack_int;
}


use std::vec;

use common::*;
use vector_monomial::VectorMonomial;
use monomial::Monomial;
use polynomial::PolyWithBorrowedMons;
use mesh::{Mesh, OShape, SideFace};
use dense_matrix::DenseMatrix;

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

pub struct WeakGradientSolver<M> {

  basis_vmons: ~[VectorMonomial<M>], // ordered ascending by monomial component dimension then monomial exponents tuple

  wgrad_comp_mons: ~[M],

  ips_basis_vmon_vs_basis_vmon_by_oshape: ~[DenseMatrix]

}

fn solve_linear_sys(matrix: &DenseMatrix, rhs: &[R]) -> ~[R] {
  // TODO: call blas routine here
  ~[0 as R]
}

impl <M:Monomial> WeakGradientSolver<M> {

  pub fn new<MESHT:Mesh<M>>(wgrad_mons_max_deg: Deg, mesh: &MESHT) -> WeakGradientSolver<M> {
    let vmons: ~[VectorMonomial<M>] = VectorMonomial::vector_mons_of_deg_le(wgrad_mons_max_deg);
    let comp_mons: ~[M] = Monomial::mons_of_deg_le(wgrad_mons_max_deg);
    let ips = make_vmon_ips_by_oshape(vmons, mesh);
    WeakGradientSolver { 
      basis_vmons: vmons,
      wgrad_comp_mons: comp_mons,
      ips_basis_vmon_vs_basis_vmon_by_oshape: ips
    }
  }
  
 /*
  * Obtain the weak gradient on the interior of the indicated oriented shape as a vector of interior-relative polynomials,
  * for a weak function on the mesh which is a monomial on the interior and 0 elsewhere.
  */
  pub fn wgrad_int_mon<MESHT:Mesh<M>>(&self, mon: M, oshape: OShape, mesh: &MESHT) -> ~[PolyWithBorrowedMons<M>] {
    let sys_rhs = vec::from_fn(self.basis_vmons.len(), |i| {
      self.wgrad_def_rhs_for_int_mon(mon.clone(), oshape, &self.basis_vmons[i], mesh)
    });
    let sol_coefs = solve_linear_sys(&self.ips_basis_vmon_vs_basis_vmon_by_oshape[*oshape], sys_rhs);
    sol_coefs.chunk_iter(self.wgrad_comp_mons.len())
             .map(|comp_coefs| PolyWithBorrowedMons { coefs: comp_coefs.to_owned(), mons: self.wgrad_comp_mons })
             .collect()
  }

 /*
  * Obtain the weak gradient on the interior of the indicated oriented shape as a vector of interior-relative polynomials,
  * for a weak function on the mesh which is a monomial on the given side and 0 elsewhere.
  */
  pub fn wgrad_side_mon<MESHT:Mesh<M>>(&self, mon: M, oshape: OShape, side_face: SideFace, mesh: &MESHT) -> ~[PolyWithBorrowedMons<M>] {
    let sys_rhs = vec::from_fn(self.basis_vmons.len(), |i| {
      self.wgrad_def_rhs_for_side_mon(mon.clone(), oshape, side_face, &self.basis_vmons[i], mesh)
    });
    let sol_coefs = solve_linear_sys(&self.ips_basis_vmon_vs_basis_vmon_by_oshape[*oshape], sys_rhs);
    sol_coefs.chunk_iter(self.wgrad_comp_mons.len())
             .map(|comp_coefs| PolyWithBorrowedMons { coefs: comp_coefs.to_owned(), mons: self.wgrad_comp_mons })
             .collect()
  }

 /*
  * Compute one component of the right hand side of the equation (WGRAD_DEF),
  *          -(v_0, div q)_T + <v_b, q.n>_bnd(T),
  * on a reference finite element of the mesh for weak function v. Here v
  * is specified as a monomial or polynomial and the supporting face on which
  * it takes the monomial or polynomial value.
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

}

fn make_vmon_ips_by_oshape<M:Monomial, MESH:Mesh<M>>(vmons: &[VectorMonomial<M>], mesh: &MESH) -> ~[DenseMatrix] {
  vec::from_fn(mesh.num_oriented_element_shapes(), |os| {
    DenseMatrix::symmetric_from_lower_triangular_fn(vmons.len(), vmons.len(), |i,j| {
      if vmons[i].mon_dim != vmons[j].mon_dim { 0 as R }
      else {
        mesh.intg_facerel_mon_on_oshape_int(vmons[i].mon * vmons[j].mon, OShape(os as u32))
      }
    })
  })
}


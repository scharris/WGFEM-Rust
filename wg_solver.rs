use common::{R};
use monomial::Monomial;
use polynomial::{Polynomial, PolyBorrowingMons};
use dense_matrix::DenseMatrix;
use mesh::{Mesh, FENum, OShape, SideFace};
use wg_basis::{WgBasis, FaceMonNum, BasisElNum};
use projection::Projector;
use variational_bilinear_form::VariationalBilinearForm;
use lapack;

use std::hashmap::HashMap;
use vbf_laplace::VBFLaplace; // TODO: remove when Rust ICE is fixed


/* METHOD
 * Let {b_i}_i be a basis for V_h^0(Omega), and vbf the bilinear form for
 * the variational problem.  Then the WG approximate solution u_h satisfies
 *   vbf(u_h, v) = (f, v_0) for all v in V_h^0(Omega)
 * Since the {b_i}_i are a basis for V_H^0, this holds iff
 *   vbf(u_h, b_i) = (f, (b_i)_0) for all basis elements b_i
 *
 * With
 *   u_h = sum_j{eta_j b_j} + Q_b g, where Q_b is L2 projection on each segment of
 * the outside boundary of Omega and 0 elsewhere, this becomes
 *   vbf(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i,
 * ie.,
 *   (sys)
 *         sum_j{ vbf(b_j, b_i) eta_j } = (f, (b_i)_0) - vbf(Q_b g, b_i) for all i
 *
 * which is a linear system we can solve for the unknown eta_j coefficients defining u_h.
 * Note that the matrix for the system m is given by m_{i,j} = vbf(b_j, b_i), which is
 * the *transpose* of the basis-vs-basis vbf values matrix.
 */

struct WgSolution<'self,Mon,MeshT> {
  basis_coefs: ~[R],
  basis: &'self WgBasis<Mon,MeshT>,
  bnd_projs: BoundaryProjections<'self,Mon>
}

/*
pub fn solve<'self, Mon: Monomial, MeshT: Mesh<Mon>, VBF: VariationalBilinearForm<'self, Mon, MeshT>>
       (vbf: &'self VBF, f: &fn(&[R])->R, g: &fn(&[R])->R) -> WgSolution<'self,Mon,MeshT> {
*/
fn solve_laplace<'a, Mon:Monomial, MeshT:Mesh<Mon>>
   (basis: &'a WgBasis<Mon,MeshT>, f: &fn(&[R])->R, g: &fn(&[R])->R) -> WgSolution<'a,Mon,MeshT> {
  
  // TODO: Replace with vbf trait parameter when Rust internal compiler error is fixed. https://github.com/mozilla/rust/issues/10201
  let vbf = &VBFLaplace::new(None, basis);

  let bnd_projs = boundary_projections(g, basis);

  let sys_m = vbf.basis_els_vs_basis_els_transpose();

  let sys_rhs = DenseMatrix::from_fn(basis.num_els(), 1, |i,_| 
    ip_on_ints(|x|f(x), BasisElNum(i), basis)
    - 
    vbf_bnd_projs_vs_bel(vbf, &bnd_projs, BasisElNum(i), basis)
  );

  let sol = lapack::solve_sparse_structurally_symmetric(&sys_m, &sys_rhs, vbf.is_symmetric());

  WgSolution {
    basis_coefs: sol,
    basis: basis,
    bnd_projs: bnd_projs,
  }
}

fn ip_on_ints<Mon:Monomial, MeshT: Mesh<Mon>>
   (f: &fn(&[R])->R,
    bel: BasisElNum,
    basis: &WgBasis<Mon,MeshT>) -> R {
  if basis.is_int_supported(bel) { 0 as R }
  else {
    let (bel_fe, bel_mon) = (basis.support_int_fe_num(bel), basis.int_mon(bel));
    basis.mesh().intg_global_fn_x_facerel_mon_on_fe_int(f, bel_mon, bel_fe)
  }
}

// Evaluate the variational bilinear form for the projection of the boundary value function g
// onto outside boundary segments vs. the given basis element. This is the vbf(Q_b g, b_i) term
// of the right hand side of (sys). The implementation uses the Element Summability and Locality
// properties of supported variational forms (see VBF module for a discussion).
fn vbf_bnd_projs_vs_bel<'a, Mon:Monomial, MeshT: Mesh<Mon>>
   (vbf: &VBFLaplace<'a,Mon,MeshT>,
    bnd_projs: &BoundaryProjections<'a,Mon>,
    bel: BasisElNum,
    basis: &WgBasis<Mon,MeshT>) -> R {

  let mesh = basis.mesh();
  
  if basis.is_int_supported(bel) {
    // Only boundary sides which are included in the bel's supporting fe can contribute.
    let fe = basis.support_int_fe_num(bel);
    let fe_oshape = mesh.oriented_shape_for_fe(fe);
    let bel_monn = basis.int_rel_mon_num(bel);
    // Sum of vbf(p,bel) for all boundary side projections p on the boundary sides of bel's fe.
    range(0, mesh.num_side_faces_for_shape(fe_oshape)).fold(0 as R, |sum, sf| {
        match bnd_projs.find(&(fe, SideFace(sf))) {
          Some(bnd_proj) => { 
            bnd_proj.foldl_numbered_terms(sum, |sum, (proj_monn, proj_coef, _)| {
              sum + proj_coef * vbf.side_mon_vs_int_mon(fe_oshape, 
                                                        FaceMonNum(proj_monn), SideFace(sf), // bnd proj mon and sf
                                                        bel_monn)                            // int basis el
            })
          }
          None => sum
        }
      })
  }
  else { // Side supported - only projections on boundary sides on one of the two including fe's can contribute.
    let bel_monn = basis.side_rel_mon_num(bel);
    let fe_incls_of_bel_supp = {
      let incls = basis.fe_inclusions_of_side_support(bel);
      [(incls.fe1, incls.side_face_in_fe1),
       (incls.fe2, incls.side_face_in_fe2)]
    };
    // Sum contributions from projections on boundary sides (if any) of both including fe's against the basis element.
    fe_incls_of_bel_supp.iter().fold(0 as R, |sum, &(fe, bel_sf)| {
      let fe_oshape = mesh.oriented_shape_for_fe(fe);
      range(0, mesh.num_side_faces_for_shape(fe_oshape)).fold(sum, |sum, sf| {
        match bnd_projs.find(&(fe, SideFace(sf))) {
          Some(bnd_proj) =>
            bnd_proj.foldl_numbered_terms(sum, |sum, (proj_monn, proj_coef, _)| {
              sum + proj_coef * vbf.side_mon_vs_side_mon_fe_contr(fe_oshape, 
                                                                  FaceMonNum(proj_monn), SideFace(sf), // bnd proj mon and side face
                                                                  bel_monn, bel_sf)                    // basis el mon and side face
            }),
          None => sum
        }
      })
    })
  }
}

type BoundaryProjections<'self,Mon> = HashMap<(FENum,SideFace), PolyBorrowingMons<'self,Mon>>;

fn boundary_projections<'a, Mon: Monomial, MeshT: Mesh<Mon>>
   (g: &fn(&[R])->R, basis: &'a WgBasis<Mon,MeshT>) -> BoundaryProjections<'a, Mon> {
 
  let mut projector = Projector::new(basis);

  let mut projs_by_fe_sf = HashMap::with_capacity(basis.mesh().num_boundary_sides());
  
  let b_fes_by_oshape_sf = basis.mesh().boundary_fes_by_oshape_side();

  for (fes_oshape, oshape_b_fes_by_sf) in b_fes_by_oshape_sf.iter().enumerate() {
    for (sf, fes) in oshape_b_fes_by_sf.iter().enumerate() {
      let fes_projs = projector.projs_to_span_fes_side_supp_basis_els(|x|g(x), *fes, OShape(fes_oshape), SideFace(sf));
      for (proj_poly, fe) in fes_projs.move_iter().zip(fes.iter()) {
        projs_by_fe_sf.insert((*fe, SideFace(sf)), proj_poly);
      }
    }
  }
  
  projs_by_fe_sf
}



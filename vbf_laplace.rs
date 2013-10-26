use common::{R};
use monomial::Monomial;
use polynomial::{PolyBorrowingMons};
use dense_matrix::DenseMatrix;
use mesh::{Mesh, OShape, SideFace};
use wg_basis::{WgBasis, FaceMonNum};
use weak_gradient::{WeakGrad, WeakGradOps};
use projection::Projector;
use variational_bilinear_form::VariationalBilinearForm;

use std::vec;
use std::cast;
use std::iter::{AdditiveIterator};

/*  Definition
 *  a_s(v, w) = sum_T{ (a wgrad_T v, wgrad_T w)_T } + sum_T{ (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T) }
 *    for piecewise functions v and w on the mesh, where v_T0 and w_T0 are the extensions of the interior
 *    polynomials v|T0 and w|T0 on T0 to the boundary of T, and Q_b is per-side projection onto the span
 *    of basis elements supported on the boundary of T.
 */

struct LaplaceVBF<'self,Mon,MeshT> {

  left_wgrad_multiplier: Option<DenseMatrix>, // post-multiplier matrix for left weak gradient in inner product

  int_mon_side_projs: ~[~[~[PolyBorrowingMons<'self,Mon>]]], // indexed by fe oshape, side face, int mon num

  weak_grad_ops: WeakGradOps<Mon>,
  
  basis: &'self WgBasis<Mon,MeshT>,

}

impl<'self, Mon:Monomial, MeshT:Mesh<Mon>> LaplaceVBF<'self,Mon,MeshT> {

  pub fn new(left_wgrad_multiplier: Option<DenseMatrix>,
             projector: &'self mut Projector<'self,Mon,MeshT>) -> LaplaceVBF<'self,Mon,MeshT> {
    
    let basis = projector.basis();

    let int_mon_side_projs = vec::from_fn(basis.mesh().num_oriented_element_shapes(), |os| {
      vec::from_fn(basis.mesh().num_side_faces_for_shape(OShape(os)), |sf| {
        projector.proj_int_mons_to_span_oshape_side_supp_basis_els(basis.ref_int_mons(), OShape(os), SideFace(sf))
      })
    });

    LaplaceVBF { 
      left_wgrad_multiplier: left_wgrad_multiplier,
      basis: basis,
      int_mon_side_projs: int_mon_side_projs,
      weak_grad_ops: basis.new_weak_grad_ops()
    }
  }

  #[inline]
  fn ip_wgrads_term(&self, wgrad_1: &WeakGrad, wgrad_2: &WeakGrad, oshape: OShape) -> R {
    let wgrads_prod = unsafe {
      match self.left_wgrad_multiplier {
        Some(ref m) => cast::transmute_mut(self).weak_grad_ops.mdot(m, wgrad_1, wgrad_2),
        None => cast::transmute_mut(self).weak_grad_ops.dot(wgrad_1, wgrad_2)
      }
    };
    self.basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, oshape)
  }

}

impl<'self, Mon:Monomial, MeshT:Mesh<Mon>> VariationalBilinearForm<'self,Mon,MeshT>
                                       for LaplaceVBF<'self,Mon,MeshT> {

  #[inline]
  fn basis(&self) -> &'self WgBasis<Mon,MeshT> {
    self.basis
  }

  #[inline]
  fn is_symmetric(&self) -> bool {
    self.left_wgrad_multiplier.is_none()
  }

  #[inline]
  fn int_mon_vs_int_mon(&self,
                        oshape: OShape,
                        monn_1: FaceMonNum,
                        monn_2: FaceMonNum) -> R {
    let (basis, mesh) = (self.basis(), self.basis().mesh());

    let ip_wgrads_term = self.ip_wgrads_term(basis.int_mon_wgrad(monn_1, oshape),
                                             basis.int_mon_wgrad(monn_2, oshape),
                                             oshape);

    // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
    // The values of v and w on the boundary are 0, leaving only
    // s(v,w) = (1/h_T) <Q_b v_T0, Q_b w_T0>_bnd(T)
    //        = (1/h_T) sum_{s in sides(T)} {<Q_b v_T0, Q_b w_T0>_s}
    let h_inv = mesh.shape_diameter_inv(oshape);
    let stab_term = h_inv * range(0, mesh.num_side_faces_for_shape(oshape)).map(|sf| {
      let projs = &self.int_mon_side_projs[*oshape][sf];
      mesh.intg_facerel_poly_x_facerel_poly_on_oshape_side(&projs[*monn_1], &projs[*monn_2], oshape, SideFace(sf))
    }).sum();

    ip_wgrads_term + stab_term
  }

  #[inline]
  fn side_mon_vs_int_mon(&self,
                         oshape: OShape,
                         side_monn: FaceMonNum, side_face: SideFace,
                         int_monn: FaceMonNum) -> R {
    let (basis, mesh) = (self.basis(), self.basis().mesh());
    
    let ip_wgrads_term = self.ip_wgrads_term(basis.side_mon_wgrad(side_monn, oshape, side_face),
                                             basis.int_mon_wgrad(int_monn, oshape),
                                             oshape);

    // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
    // With v being our side supported basis element and w the interior supported element,
    // s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
    //        = (1/h_T) <-v, Q_b w_T0>_bnd(T)
    //        = (1/h_T) <-v, Q_b w_T0>_s where s is the supporting side face of v
    let stab_term = {
      let side_mon = basis.side_mons_for_oshape_side(oshape, side_face)[*side_monn].clone();
      let int_proj = &self.int_mon_side_projs[*oshape][*side_face][*int_monn];
      let ip = -mesh.intg_facerel_mon_x_facerel_poly_on_oshape_side(side_mon, int_proj, oshape, side_face);
      let h_inv = mesh.shape_diameter_inv(oshape);
      h_inv * ip
    };

    ip_wgrads_term + stab_term
  }
  
  #[inline]
  fn int_mon_vs_side_mon(&self,
                         oshape: OShape,
                         int_monn: FaceMonNum,
                         side_monn: FaceMonNum, side_face: SideFace) -> R {
    let (basis, mesh) = (self.basis(), self.basis().mesh());

    let ip_wgrads_term = self.ip_wgrads_term(basis.int_mon_wgrad(int_monn, oshape),
                                             basis.side_mon_wgrad(side_monn, oshape, side_face),
                                             oshape);

    // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
    // With v being our interior supported basis element and w the side supported element, 
    // s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
    //        = (1/h_T) <Q_b v_T0, -w>_bnd(T)
    //        = (1/h_T) <Q_b v_T0, -w>_s where s is the supporting side face of w
    let stab_term = {
      let side_mon = basis.side_mons_for_oshape_side(oshape, side_face)[*side_monn].clone();
      let int_proj = &self.int_mon_side_projs[*oshape][*side_face][*int_monn];
      let ip = -mesh.intg_facerel_mon_x_facerel_poly_on_oshape_side(side_mon, int_proj, oshape, side_face);
      let h_inv = mesh.shape_diameter_inv(oshape);
      h_inv * ip
    };

    ip_wgrads_term + stab_term
  }

  #[inline]
  fn side_mon_vs_side_mon_fe_contr(&self,
                                   oshape: OShape,
                                   monn_1: FaceMonNum, side_face_1: SideFace,
                                   monn_2: FaceMonNum, side_face_2: SideFace) -> R {
    let (basis, mesh) = (self.basis(), self.basis().mesh());
    
    let ip_wgrads_term = self.ip_wgrads_term(basis.side_mon_wgrad(monn_1, oshape, side_face_1),
                                             basis.side_mon_wgrad(monn_2, oshape, side_face_2),
                                             oshape);

    // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
    // With v and w being our side supported elements, the projections of the interior extensions are 0, so
    // s(v,w) = (1/h_T) <-v, -w>_bnd(T)
    //        = (1/h_T) sum_{s in sides(T)} <v, w>_s
    //        = | (1/h_T) <v, w>_s,  if v and w have a common support side face s
    //          | 0, otherwise
    let stab_term =
      if side_face_1 != side_face_2 { 0 as R }
      else {
        let common_supp_side = side_face_1;
        let side_mons = basis.side_mons_for_oshape_side(oshape, common_supp_side);
        let ip = mesh.intg_facerel_mon_on_oshape_side(side_mons[*monn_1] * side_mons[*monn_2], oshape, common_supp_side);
        mesh.shape_diameter_inv(oshape) * ip
      };

    ip_wgrads_term + stab_term
  }
}


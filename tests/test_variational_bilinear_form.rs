use variational_bilinear_form::VariationalBilinearForm;
use common::{R};
use monomial::{Monomial, Mon2d, MaxMonDeg};
use mesh::{Mesh, OShape, SideFace};
use rectangle_mesh::{RectMesh, MeshCoord};
use wg_basis::{WGBasis, FaceMonNum, BasisElNum};


struct AsymmetricTestVBF<'self,Mon,MeshT> {
  basis: &'self WGBasis<Mon,MeshT>,
}

impl<'self,Mon:Monomial,MeshT:Mesh<Mon>> VariationalBilinearForm<'self,Mon,MeshT>
                                     for AsymmetricTestVBF<'self,Mon,MeshT> {

  fn basis(&self) -> &'self WGBasis<Mon,MeshT> {
    self.basis
  }

  fn is_symmetric(&self) -> bool { false }

  fn int_mon_vs_int_mon(&self,
                        oshape: OShape,
                        monn_1: FaceMonNum,
                        monn_2: FaceMonNum) -> R {
    (&[*oshape, *monn_1, *monn_2]).hash() as R
  }

  fn side_mon_vs_int_mon(&self,
                         oshape: OShape,
                         side_monn: FaceMonNum, side_face: SideFace,
                         int_monn: FaceMonNum) -> R {
    (&[*oshape, *side_monn, *side_face, *int_monn]).hash() as R
  }
  
  fn int_mon_vs_side_mon(&self,
                         oshape: OShape,
                         int_monn: FaceMonNum,
                         side_monn: FaceMonNum, side_face: SideFace) -> R {
    (&[*oshape, *int_monn, *side_monn, *side_face]).hash() as R
  }

  fn side_mon_vs_side_mon_fe_contr(&self,
                                   oshape: OShape,
                                   monn_1: FaceMonNum, side_face_1: SideFace,
                                   monn_2: FaceMonNum, side_face_2: SideFace) -> R {
    (&[*oshape, *monn_1, *side_face_1, *monn_2, *side_face_2]).hash() as R
  }
}


struct SymmetricTestVBF<'self, Mon, MeshT> {
  basis: &'self WGBasis<Mon,MeshT>,
}

impl<'self,Mon:Monomial,MeshT:Mesh<Mon>> VariationalBilinearForm<'self,Mon,MeshT>
                                      for SymmetricTestVBF<'self,Mon,MeshT> {

  fn basis(&self) -> &'self WGBasis<Mon,MeshT> {
    self.basis
  }

  fn is_symmetric(&self) -> bool { true }

  fn int_mon_vs_int_mon(&self,
                        oshape: OShape,
                        monn_1: FaceMonNum,
                        monn_2: FaceMonNum) -> R {
    (&[*oshape, *monn_1, *monn_2]).hash() as R +
    (&[*oshape, *monn_2, *monn_1]).hash() as R
  }

  fn side_mon_vs_int_mon(&self,
                         oshape: OShape,
                         side_monn: FaceMonNum, side_face: SideFace,
                         int_monn: FaceMonNum) -> R {
    (&[*oshape, *side_monn, *side_face, *int_monn]).hash() as R +
    (&[*oshape, *int_monn, *side_monn, *side_face]).hash() as R
  }
  
  fn int_mon_vs_side_mon(&self,
                         oshape: OShape,
                         int_monn: FaceMonNum,
                         side_monn: FaceMonNum, side_face: SideFace) -> R {
    (&[*oshape, *side_monn, *side_face, *int_monn]).hash() as R +
    (&[*oshape, *int_monn, *side_monn, *side_face]).hash() as R
  }

  fn side_mon_vs_side_mon_fe_contr(&self,
                                   oshape: OShape,
                                   monn_1: FaceMonNum, side_face_1: SideFace,
                                   monn_2: FaceMonNum, side_face_2: SideFace) -> R {
    (&[*oshape, *monn_1, *side_face_1, *monn_2, *side_face_2]).hash() as R +
    (&[*oshape, *monn_2, *side_face_2, *monn_1, *side_face_1]).hash() as R
  }
}

#[test]
fn test_asymmetric() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[2.,2.], ~[MeshCoord(5),MeshCoord(4)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));

  let vbf: AsymmetricTestVBF<Mon2d,RectMesh<Mon2d>> = AsymmetricTestVBF { basis: basis };

  let m = vbf.basis_els_vs_basis_els_transpose();

  for r in range(0, basis.num_els()) {
    for c in range(0, basis.num_els()) {
      let (r,c) = (BasisElNum(r), BasisElNum(c));

      if basis.is_int_supported(r) {
        let r_int_fe = basis.support_int_fe_num(r);
        let r_monn = basis.int_rel_mon_num(r);
        if basis.is_int_supported(c) { // (r,c) is an int vs int el pair
          let c_int_fe = basis.support_int_fe_num(c);
          let c_monn = basis.int_rel_mon_num(c);
          let ip = if r_int_fe == c_int_fe { vbf.int_mon_vs_int_mon(OShape(0), c_monn, r_monn) } else { 0 as R };
          assert_eq!(m.get(*r,*c), ip);
        }
        else { // (r,c) is an int vs side el pair
          let c_incls = basis.fe_inclusions_of_side_support(c);
          let c_monn = basis.side_rel_mon_num(c);
          let ip = if c_incls.fe1 == r_int_fe || c_incls.fe2 == r_int_fe {
            let c_sf = if c_incls.fe1 == r_int_fe { c_incls.side_face_in_fe1 } else { c_incls.side_face_in_fe2 };
            vbf.side_mon_vs_int_mon(OShape(0), c_monn, c_sf, r_monn)
          } else { 0 as R };
          assert_eq!(m.get(*r,*c), ip);
        }
      }
      else { // first basis element is side-supported
        let r_incls = basis.fe_inclusions_of_side_support(r);
        let r_monn = basis.side_rel_mon_num(r);

        if basis.is_int_supported(c) { // (r,c) is an side vs int el pair
          let c_int_fe = basis.support_int_fe_num(c);
          let c_monn = basis.int_rel_mon_num(c);
          let ip = if r_incls.fe1 == c_int_fe || r_incls.fe2 == c_int_fe {
            let r_sf = if r_incls.fe1 == c_int_fe { r_incls.side_face_in_fe1 } else { r_incls.side_face_in_fe2 };
            vbf.int_mon_vs_side_mon(OShape(0), c_monn, r_monn, r_sf)
          } else { 0 as R };
          assert_eq!(m.get(*r,*c), ip);
        }
        else { // (r,c) is a side vs side el pair
          let c_incls = basis.fe_inclusions_of_side_support(c);
          let c_monn = basis.side_rel_mon_num(c);

          let mut ip = 0 as R;
          
          if r_incls.fe1 == c_incls.fe1 {
            ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe1, r_monn, r_incls.side_face_in_fe1);
          } 
          else if r_incls.fe1 == c_incls.fe2 {
            ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe2, r_monn, r_incls.side_face_in_fe1);
          }
          
          if r_incls.fe2 == c_incls.fe1 {
            ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe1, r_monn, r_incls.side_face_in_fe2);
          } 
          else if r_incls.fe2 == c_incls.fe2 {
            ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe2, r_monn, r_incls.side_face_in_fe2);
          }
          assert_eq!(m.get(*r, *c), ip);
        }
      }
    }
  }
}

#[test]
fn test_symmetric() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[2.,2.], ~[MeshCoord(5),MeshCoord(4)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));

  let vbf: SymmetricTestVBF<Mon2d,RectMesh<Mon2d>> = SymmetricTestVBF { basis: basis };

  let m = vbf.basis_els_vs_basis_els_transpose();

  for r in range(0, basis.num_els()) {
    for c in range(0, basis.num_els()) {
      let (r,c) = (BasisElNum(r), BasisElNum(c));

      if basis.is_int_supported(r) {
        let r_int_fe = basis.support_int_fe_num(r);
        let r_monn = basis.int_rel_mon_num(r);
        if basis.is_int_supported(c) { // (r,c) is an int vs int el pair
          let c_int_fe = basis.support_int_fe_num(c);
          let c_monn = basis.int_rel_mon_num(c);
          let ip = if r <= c && r_int_fe == c_int_fe { vbf.int_mon_vs_int_mon(OShape(0), c_monn, r_monn) } else { 0 as R };
          assert_eq!(m.get(*r,*c), ip);
        }
        else { // (r,c) is an int vs side el pair
          let c_incls = basis.fe_inclusions_of_side_support(c);
          let c_monn = basis.side_rel_mon_num(c);
          let ip = if r <= c && (c_incls.fe1 == r_int_fe || c_incls.fe2 == r_int_fe) {
            let c_sf = if c_incls.fe1 == r_int_fe { c_incls.side_face_in_fe1 } else { c_incls.side_face_in_fe2 };
            vbf.side_mon_vs_int_mon(OShape(0), c_monn, c_sf, r_monn)
          } else { 0 as R };
          assert_eq!(m.get(*r,*c), ip);
        }
      }
      else { // first basis element is side-supported
        let r_incls = basis.fe_inclusions_of_side_support(r);
        let r_monn = basis.side_rel_mon_num(r);

        if basis.is_int_supported(c) { // (r,c) is an side vs int el pair
          let c_int_fe = basis.support_int_fe_num(c);
          let c_monn = basis.int_rel_mon_num(c);
          let ip = if r <= c && (r_incls.fe1 == c_int_fe || r_incls.fe2 == c_int_fe) {
            let r_sf = if r_incls.fe1 == c_int_fe { r_incls.side_face_in_fe1 } else { r_incls.side_face_in_fe2 };
            vbf.int_mon_vs_side_mon(OShape(0), c_monn, r_monn, r_sf)
          } else { 0 as R };
          assert_eq!(m.get(*r,*c), ip);
        }
        else { // (r,c) is a side vs side el pair
          let c_incls = basis.fe_inclusions_of_side_support(c);
          let c_monn = basis.side_rel_mon_num(c);
          let mut ip = 0 as R;
          if r <= c {
            if r_incls.fe1 == c_incls.fe1 {
              ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe1, r_monn, r_incls.side_face_in_fe1);
            } 
            else if r_incls.fe1 == c_incls.fe2 {
              ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe2, r_monn, r_incls.side_face_in_fe1);
            }
            if r_incls.fe2 == c_incls.fe1 {
              ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe1, r_monn, r_incls.side_face_in_fe2);
            } 
            else if r_incls.fe2 == c_incls.fe2 {
              ip += vbf.side_mon_vs_side_mon_fe_contr(OShape(0), c_monn, c_incls.side_face_in_fe2, r_monn, r_incls.side_face_in_fe2);
            }
          }
          assert_eq!(m.get(*r, *c), ip);
        }
      }
    }
  }
}


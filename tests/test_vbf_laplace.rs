use vbf_laplace::LaplaceVBF;
use variational_bilinear_form::VariationalBilinearForm;
use common::{Deg};
use monomial::{Mon2d, MaxMonDeg};
use dense_matrix::DenseMatrix;
use mesh::{Mesh, OShape, SideFace};
use rectangle_mesh::{RectMesh, MeshCoord};
use wg_basis::{WgBasis, FaceMonNum};

use std::num::{sqrt};

#[test]
fn test_is_symmetric() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
 
  let m = DenseMatrix::from_fn(2,2, |r,c| { if r != c { 1. } else { 0. } });
  let vbf_asym = LaplaceVBF::new(Some(m), basis);
  assert!(!vbf_asym.is_symmetric());

  let vbf_sym = LaplaceVBF::new(None, basis);
  assert!(vbf_sym.is_symmetric());
}

#[test]
fn test_int_vs_int_sym() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  
  // Our interior mons will be xy and y^2, which have face numbers of 5 and 2 in the
  // reference sequence of interior monomials: 1, y, y^2, y^3, x, xy, xy^2, x^2, x^2y, x^3.
  let wgrad_int_xy = basis.int_mon_wgrad(FaceMonNum(5), OShape(0));
  let wgrad_int_y2 = basis.int_mon_wgrad(FaceMonNum(2), OShape(0));
 
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.dot(wgrad_int_xy, wgrad_int_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));
  

  // stabilization term
  // Since our interior monomials v and w are 0 on the boundary, and they project to themselves
  // on each side's space of basis elements supported on the side, we have 
  // s(v,w) = (1/h_T) sum_{s in sides(T)} {<Q_b v_T0, Q_b w_T0>_s}
  //        = (1/h_T) sum_{s in sides(T)} {<v_T0, w_T0>_s}
  let stab_term = {
    let x = Mon2d { exps: [Deg(1), Deg(0)] };
    let y = Mon2d { exps: [Deg(0), Deg(1)] };
    let sides_ip = basis.mesh().intg_facerel_mon_on_oshape_side(y*y*y, OShape(0), SideFace(1)) // right side
                 + basis.mesh().intg_facerel_mon_on_oshape_side(x, OShape(0), SideFace(3));    // top side
    let h_inv = basis.mesh().shape_diameter_inv(OShape(0));
    h_inv * sides_ip
  };
  
  let vbf = LaplaceVBF::new(None, basis);
  
  assert_eq!(vbf.int_mon_vs_int_mon(OShape(0), FaceMonNum(5), FaceMonNum(2)),
             ips_term + stab_term);
}

#[test]
fn test_int_vs_int_asym() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  
  // Our interior mons will be xy and y^2, which have face numbers of 5 and 2 in the
  // reference sequence of interior monomials: 1, y, y^2, y^3, x, xy, xy^2, x^2, x^2y, x^3.
  let wgrad_int_xy = basis.int_mon_wgrad(FaceMonNum(5), OShape(0));
  let wgrad_int_y2 = basis.int_mon_wgrad(FaceMonNum(2), OShape(0));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_int_xy, wgrad_int_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term
  // Since our interior monomials v and w are 0 on the boundary, and they project to themselves
  // on each side's space of basis elements supported on the side, we have 
  // s(v,w) = (1/h_T) sum_{s in sides(T)} {<Q_b v_T0, Q_b w_T0>_s}
  //        = (1/h_T) sum_{s in sides(T)} {<v_T0, w_T0>_s}
  let stab_term = {
    let x = Mon2d { exps: [Deg(1), Deg(0)] };
    let y = Mon2d { exps: [Deg(0), Deg(1)] };
    let sides_ip = basis.mesh().intg_facerel_mon_on_oshape_side(y*y*y, OShape(0), SideFace(1)) // right side
                 + basis.mesh().intg_facerel_mon_on_oshape_side(x, OShape(0), SideFace(3));    // top side
    let h_inv = basis.mesh().shape_diameter_inv(OShape(0));
    h_inv * sides_ip
  };
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.int_mon_vs_int_mon(OShape(0), FaceMonNum(5), FaceMonNum(2)),
             ips_term + stab_term);
}

#[test]
fn test_top_side_vs_int() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be x on the top side and xy on the interior.
  //   Reference side monomials for horizontal sides: 1, x, x^2.
  //   Reference interior monomials: 1, y, y^2, y^3, x, xy, xy^2, x^2, x^2y, x^3.
  let wgrad_top_x = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(3));
  let wgrad_int_xy = basis.int_mon_wgrad(FaceMonNum(5), OShape(0));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_top_x, wgrad_int_xy);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v being our side supported basis element and w the interior supported element,
  // s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  //        = (1/h_T) <-v, Q_b w_T0>_bnd(T)
  //        = (1/h_T) <-v, Q_b w_T0>_s where s is the supporting side face of v.
  let stab_term = {
    let x = Mon2d { exps: [Deg(1), Deg(0)] };
    let sides_ip = -basis.mesh().intg_facerel_mon_on_oshape_side(x*x, OShape(0), SideFace(3)); // top side
    let h_inv = basis.mesh().shape_diameter_inv(OShape(0));
    h_inv * sides_ip
  };
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.side_mon_vs_int_mon(OShape(0), FaceMonNum(1), SideFace(3), FaceMonNum(5)),
             ips_term + stab_term);
}

#[test]
fn test_bottom_side_vs_int() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be x on the bottom side and xy on the interior.
  //   Reference side monomials for horizontal sides: 1, x, x^2.
  //   Reference interior monomials: 1, y, y^2, y^3, x, xy, xy^2, x^2, x^2y, x^3.
  let wgrad_bottom_x = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(2));
  let wgrad_int_xy = basis.int_mon_wgrad(FaceMonNum(5), OShape(0));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_bottom_x, wgrad_int_xy);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v being our side supported basis element and w the interior supported element,
  // s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  //        = (1/h_T) <-v, Q_b w_T0>_bnd(T)
  //        = (1/h_T) <-v, Q_b w_T0>_s where s is the supporting side face of v.
  let stab_term = 0.; // Because the xy projection is 0 on the right side of the inner product.
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.side_mon_vs_int_mon(OShape(0), FaceMonNum(1), SideFace(2), FaceMonNum(5)),
             ips_term + stab_term);
}

#[test]
fn test_int_vs_right_side() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be xy on the interior and y^2 on the right side.
  //   Reference interior monomials: 1, y, y^2, y^3, x, xy, xy^2, x^2, x^2y, x^3.
  //   Reference side monomials for vertical sides: 1, y, y^2.
  let wgrad_int_xy = basis.int_mon_wgrad(FaceMonNum(5), OShape(0));
  let wgrad_right_y2 = basis.side_mon_wgrad(FaceMonNum(2), OShape(0), SideFace(1));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_int_xy, wgrad_right_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v being our interior supported basis element and w the side supported element, 
  // s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  //        = (1/h_T) <Q_b v_T0, -w>_bnd(T)
  //        = (1/h_T) <Q_b v_T0, -w>_s where s is the supporting face of w.
  let stab_term = {
    let y = Mon2d { exps: [Deg(0), Deg(1)] };
    let sides_ip = -basis.mesh().intg_facerel_mon_on_oshape_side(y*y*y, OShape(0), SideFace(1)); // right side
    let h_inv = basis.mesh().shape_diameter_inv(OShape(0));
    h_inv * sides_ip
  };
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.int_mon_vs_side_mon(OShape(0), FaceMonNum(5), FaceMonNum(2), SideFace(1)),
             ips_term + stab_term);
}

#[test]
fn test_int_vs_left_side() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be xy on the interior and y^2 on the left side.
  //   Reference interior monomials: 1, y, y^2, y^3, x, xy, xy^2, x^2, x^2y, x^3.
  //   Reference side monomials for vertical sides: 1, y, y^2.
  let wgrad_int_xy = basis.int_mon_wgrad(FaceMonNum(5), OShape(0));
  let wgrad_left_y2 = basis.side_mon_wgrad(FaceMonNum(2), OShape(0), SideFace(0));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_int_xy, wgrad_left_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v being our interior supported basis element and w the side supported element, 
  // s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  //        = (1/h_T) <Q_b v_T0, -w>_bnd(T)
  //        = (1/h_T) <Q_b v_T0, -w>_s where s is the supporting face of w.
  let stab_term = 0.; // Because the xy projection is 0 on the left side of the inner product.
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.int_mon_vs_side_mon(OShape(0), FaceMonNum(5), FaceMonNum(2), SideFace(0)),
             ips_term + stab_term);
}

#[test]
fn test_right_side_vs_left_side() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be y on the right side and y^2 on the left side.
  //   Reference side monomials for vertical sides: 1, y, y^2.
  let wgrad_right_y = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(1));
  let wgrad_left_y2 = basis.side_mon_wgrad(FaceMonNum(2), OShape(0), SideFace(0));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_right_y, wgrad_left_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v and w being our side supported elements, the projections of the interior extensions are 0, so
  // s(v,w) = (1/h_T) <-v, -w>_bnd(T)
  //        = (1/h_T) sum_{s in sides(T)} <v, w>_s
  //        = | (1/h_T) <v, w>_s,  if v and w have a common support side face s
  //          | 0, otherwise
  let stab_term = 0.; 
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.side_mon_vs_side_mon_fe_contr(OShape(0), FaceMonNum(1), SideFace(1), FaceMonNum(2), SideFace(0)),
             ips_term + stab_term);
}

#[test]
fn test_bottom_side_vs_left_side() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be x on the bottom side and y^2 on the left side.
  //   Reference side monomials for horizontal sides: 1, x, x^2.
  //   Reference side monomials for vertical sides: 1, y, y^2.
  let wgrad_bottom_x = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(2));
  let wgrad_left_y2 = basis.side_mon_wgrad(FaceMonNum(2), OShape(0), SideFace(0));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_bottom_x, wgrad_left_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v and w being our side supported elements, the projections of the interior extensions are 0, so
  // s(v,w) = (1/h_T) <-v, -w>_bnd(T)
  //        = (1/h_T) sum_{s in sides(T)} <v, w>_s
  //        = | (1/h_T) <v, w>_s,  if v and w have a common support side face s
  //          | 0, otherwise
  let stab_term = 0.; 
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.side_mon_vs_side_mon_fe_contr(OShape(0), FaceMonNum(1), SideFace(2), FaceMonNum(2), SideFace(0)),
             ips_term + stab_term);
}

#[test]
fn test_right_side_vs_right_side() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be y on the right side and y^2 on the right side.
  //   Reference side monomials for vertical sides: 1, y, y^2.
  let wgrad_right_y = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(1));
  let wgrad_right_y2 = basis.side_mon_wgrad(FaceMonNum(2), OShape(0), SideFace(1));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_right_y, wgrad_right_y2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v and w being our side supported elements, the projections of the interior extensions are 0, so
  // s(v,w) = (1/h_T) <-v, -w>_bnd(T)
  //        = (1/h_T) sum_{s in sides(T)} <v, w>_s
  //        = | (1/h_T) <v, w>_s,  if v and w have a common support side face s
  //          | 0, otherwise
  let stab_term = {
    let y = Mon2d { exps: [Deg(0), Deg(1)] };
    let ip = basis.mesh().intg_facerel_mon_on_oshape_side(y*y*y, OShape(0), SideFace(1));
    basis.mesh().shape_diameter_inv(OShape(0)) * ip
  };
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.side_mon_vs_side_mon_fe_contr(OShape(0), FaceMonNum(1), SideFace(1), FaceMonNum(2), SideFace(1)),
             ips_term + stab_term);
}

#[test]
fn test_bottom_side_vs_bottom_side() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
 
  // Our monomials will be x on the bottom side and x^2 on the bottom side.
  //   Reference side monomials for horizontal sides: 1, x, x^2.
  let wgrad_bottom_x = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(2));
  let wgrad_bottom_x2 = basis.side_mon_wgrad(FaceMonNum(2), OShape(0), SideFace(2));

  let m = DenseMatrix::from_rows(2,2,
    [~[1./sqrt(2.), -1./sqrt(2.)],
     ~[1./sqrt(2.),  1./sqrt(2.)]]);
  
  let mut wgrad_ops = basis.new_weak_grad_ops();
  let wgrads_prod = wgrad_ops.mdot(&m, wgrad_bottom_x, wgrad_bottom_x2);
  let ips_term = basis.mesh().intg_facerel_poly_on_oshape_int(&wgrads_prod, OShape(0));

  // stabilization term: s(v,w) = (1/h_T) <Q_b v_T0 - v, Q_b w_T0 - w>_bnd(T)
  // With v and w being our side supported elements, the projections of the interior extensions are 0, so
  // s(v,w) = (1/h_T) <-v, -w>_bnd(T)
  //        = (1/h_T) sum_{s in sides(T)} <v, w>_s
  //        = | (1/h_T) <v, w>_s,  if v and w have a common support side face s
  //          | 0, otherwise
  let stab_term = {
    let x = Mon2d { exps: [Deg(1), Deg(0)] };
    let ip = basis.mesh().intg_facerel_mon_on_oshape_side(x*x*x, OShape(0), SideFace(2));
    basis.mesh().shape_diameter_inv(OShape(0)) * ip
  };
  
  let vbf = LaplaceVBF::new(Some(m), basis);
  
  assert_eq!(vbf.side_mon_vs_side_mon_fe_contr(OShape(0), FaceMonNum(1), SideFace(2), FaceMonNum(2), SideFace(2)),
             ips_term + stab_term);
}


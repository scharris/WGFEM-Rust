use projection::{Projector};
use common::{R, Deg, Dim, pow};
use polynomial::{Polynomial, PolyBorrowingMons, approx_equiv};
use monomial::{Monomial, Mon2d, MaxMonDeg}; 
use mesh::{FENum, OShape};
use rectangle_mesh::{RectMesh, MeshCoord};
use wg_basis::{WgBasis};

use std::vec;

#[test]
fn test_int_L2_inner_products_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let projector: ~Projector<Mon2d,RectMesh<Mon2d>>  = Projector::new(basis);
  
  let int_ips = projector.ips_int_mons_by_oshape[0];
  assert_eq!(int_ips.get(0,0), 1.);    // one vs one
  assert_eq!(int_ips.get(0,1), 1./2.); // one vs y
  assert_eq!(int_ips.get(0,2), 1./3.); // one vs y^2
  assert_eq!(int_ips.get(0,3), 1./2.); // one vs x
  assert_eq!(int_ips.get(0,4), 1./4.); // one vs xy
  assert_eq!(int_ips.get(0,5), 1./3.); // one vs x^2
}

#[test]
fn test_side_L2_inner_products_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let projector: ~Projector<Mon2d,RectMesh<Mon2d>>  = Projector::new(basis);

  let right_ips = projector.ips_side_mons_by_oshape_side[0][3];
  assert_eq!(right_ips.get(0,0), 1.);    // one vs one
  assert_eq!(right_ips.get(0,1), 1./2.); // one vs y
  assert_eq!(right_ips.get(1,1), 1./3.); // y vs y

  let top_ips = projector.ips_side_mons_by_oshape_side[0][1];
  assert_eq!(top_ips.get(0,0), 1.);    // one vs one
  assert_eq!(top_ips.get(0,1), 1./2.); // one vs x
  assert_eq!(top_ips.get(1,1), 1./3.); // x vs x
}


#[test]
fn test_identity_proj_fe0() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(4));
  let mut projector: ~Projector<Mon2d,RectMesh<Mon2d>>  = Projector::new(basis);

  // 0 x^0y^0 + -2 x^0y^1 + 0 x^0y^2 + -0 x^0y^3 + 1 x^1y^0 + -3 x^1y^1 + -2 x^1y^2 + 0 x^2y^0 + 1 x^2y^1 + -0 x^3y^0
  let p = PolyBorrowingMons::new(~[0., -2., 0., 0., 1., -3., -2., 0., 1., 0.], basis.int_mons());

  let g = |x: &[R]| p.value_at(x);
  
  let proj = projector.projs_to_span_int_supp_basis_els(g, &[FENum(0)], OShape(0));

  assert!(approx_equiv(&proj[0], &p, 1e-10));
}

#[test]
fn test_identity_proj_fe4() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(4)]);
 
  let fe = FENum(4);
  let int_orig_0 = rmesh.fe_interior_origin_comp(fe, Dim(0));
  let int_orig_1 = rmesh.fe_interior_origin_comp(fe, Dim(1));
  
  let basis = WgBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(4));
  let mut projector: ~Projector<Mon2d,RectMesh<Mon2d>>  = Projector::new(basis);

  let g = |x: &[R]| {
    // x^2 y - 2 x y^2
    pow(x[0] - int_orig_0, 2) * (x[1] - int_orig_1) - 2. *  pow(x[0] - int_orig_0, 1) * pow(x[1] - int_orig_1, 2)
  };
  
  let proj = projector.projs_to_span_int_supp_basis_els(g, &[fe], OShape(0));
  
  let p = PolyBorrowingMons::new(~[0., 0., 0., 0., 0., 0., -2., 0., 1., 0.], basis.int_mons()); // x^2 y - 2 x y^2

  assert!(approx_equiv(&proj[0], &p, 1e-10));
}



// TODO: Side projection tests.


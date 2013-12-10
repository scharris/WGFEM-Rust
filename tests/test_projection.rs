use projection::{Projector};
use common::{R, Deg, pow};
use polynomial::{Polynomial, PolyBorrowingMons, approx_equiv, approx_equiv_v};
use monomial::{Monomial, Mon2d, MaxMonDeg}; 
use mesh::{FENum, OShape, SideFace};
use rectangle_mesh::{RectMesh, MeshCoord};
use wg_basis::{WGBasis};
use lapack;

#[test]
fn test_do_lapack_init() {
  lapack::init(); // TODO: Do this somewhere else, where it's gauranteed to be run before other tests as part of each test setup.
}


#[test]
fn test_identity_proj_fe0() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = &WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  let mut projector: Projector<Mon2d,RectMesh<Mon2d>>  = Projector::new(basis);

  // 0 x^0y^0 + -2 x^0y^1 + 0 x^0y^2 + -0 x^0y^3 + 1 x^1y^0 + -3 x^1y^1 + -2 x^1y^2 + 0 x^2y^0 + 1 x^2y^1 + -0 x^3y^0
  let p = PolyBorrowingMons::new(~[0., -2., 0., 0., 1., -3., -2., 0., 1., 0.], basis.ref_int_mons());

  let g = |x: &[R]| p.value_at(x);
  
  let proj = projector.projs_to_span_fes_int_supp_basis_els(g, &[FENum(0)], OShape(0));

  assert!(approx_equiv(&proj[0], &p, 1e-10));
}

#[test]
fn test_int_projs() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,4.], ~[MeshCoord(3),MeshCoord(4)]);
 
  let fe4_int_orig_0 = 1.;
  let fe4_int_orig_1 = 1.;
  
  let basis = &WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  let mut projector: Projector<Mon2d,RectMesh<Mon2d>>  = Projector::with_rhs_cols_capacity(basis, 2); // force a reallocation of rhs buffer

  let g = |x: &[R]| {
    // polynomial x^2 y - 2 x y^2 relative to fe 4 interior origin
    pow(x[0] - fe4_int_orig_0, 2) * (x[1] - fe4_int_orig_1) - 2. *  pow(x[0] - fe4_int_orig_0, 1) * pow(x[1] - fe4_int_orig_1, 2)
  };
    
  let projs = projector.projs_to_span_fes_int_supp_basis_els(g, &[FENum(0),FENum(4),FENum(0),FENum(4),FENum(0)], OShape(0));
  
  // interior mons: x^0y^0, x^0y^1, x^0y^2, x^0y^3, x^1y^0, x^1y^1, x^1y^2, x^2y^0, x^2y^1, x^3y^0
  
  // g is x^2 y - 2 x y^2 expressed locally on fe4's interior.
  let fe4_proj = PolyBorrowingMons::new(~[0., 0., 0., 0., 0., 0., -2., 0., 1., 0.], basis.ref_int_mons());
  
  //  g is (x-1)^2 (y-1) - 2 (x-1)(y-1)^2 = x^2*y - 2*x*y^2 - x^2 + 2*x*y + 2*y^2 - 3*y + 1 expressed locally on fe0's interior.
  let fe0_proj = PolyBorrowingMons::new(~[1., -3., 2., 0., 0., 2., -2., -1., 1., 0.], basis.ref_int_mons());

  assert!(approx_equiv(&projs[0], &fe0_proj, 1e-10));
  assert!(approx_equiv(&projs[1], &fe4_proj, 1e-10));
  assert!(approx_equiv(&projs[2], &fe0_proj, 1e-10));
  assert!(approx_equiv(&projs[3], &fe4_proj, 1e-10));
  assert!(approx_equiv(&projs[4], &fe0_proj, 1e-10));
}

#[test]
fn test_right_side_projs() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,4.], ~[MeshCoord(3),MeshCoord(4)]);
 
  let fe5_int_orig_0 = 2.;
  let fe5_int_orig_1 = 1.;
  
  let basis = &WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  let mut projector: Projector<Mon2d,RectMesh<Mon2d>>  = Projector::with_rhs_cols_capacity(basis, 2); // force a reallocation of rhs buffer

  // The function to be projected is the polynomial whose representation relative to fe 5's interior is
  //   p_fe5(x,y) = 2x^2 - 3y^2 - 5x - 4y (relative to fe 5 origin at (2,1)).
  // which expressed relative to the global origin (or also fe 0's interior) is polynomial
  //   p(x,y) = 2(x-2)^2 - 3(y-1)^2 - 5(x-2) - 4(y-1)
  let g = |x: &[R]| {
    let (x,y) = (x[0] - fe5_int_orig_0,
                 x[1] - fe5_int_orig_1);
    2.*pow(x, 2) - 3.*pow(y, 2) - 5.*x - 4.*y
  };

  let projs = projector.projs_to_span_fes_side_supp_basis_els(g, &[FENum(0),FENum(5),FENum(0),FENum(5),FENum(0)], OShape(0), SideFace(1));

  let right_side_mons = basis.side_mons_for_fe_side(FENum(5), SideFace(1)); // vertical side monomials: x^0y^0, x^0y^1, x^0y^2

  // At a point (x,y) on fe5's right side, the polynomial to be projected has value
  //   p(x,y) = 2(1^2) - 3y^2 - 5(1) - 4y
  //          = -3 - 4y - 3y^2
  let fe5_proj = PolyBorrowingMons::new(~[-3., -4., -3.], right_side_mons);
  
  // At a point (x,y) on fe0's right side, the polynomial to be projected has value
  //   p(x,y) = 2(1-2)^2 - 3(y-1)^2 - 5(1-2) - 4(y-1)
  //          = 8 + 2y - 3y^2
  let fe0_proj = PolyBorrowingMons::new(~[8., 2., -3.], right_side_mons);

  assert!(approx_equiv(&projs[0], &fe0_proj, 1e-10));
  assert!(approx_equiv(&projs[1], &fe5_proj, 1e-10));
  assert!(approx_equiv(&projs[2], &fe0_proj, 1e-10));
  assert!(approx_equiv(&projs[3], &fe5_proj, 1e-10));
  assert!(approx_equiv(&projs[4], &fe0_proj, 1e-10));
}

#[test]
fn test_int_mons_side_projs() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[6.,12.], ~[MeshCoord(3),MeshCoord(4)]);
  let basis = &WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  let mut projector: Projector<Mon2d,RectMesh<Mon2d>>  = Projector::new(basis);

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  
  let left_side = SideFace(0);
  let right_side = SideFace(1);
  let bottom_side = SideFace(2);
  let top_side = SideFace(3);

  let vert_side_mons = basis.side_mons_for_oshape_side(OShape(0), left_side);
  let horz_side_mons = basis.side_mons_for_oshape_side(OShape(0), SideFace(2));

  // Test identity projections onto vertical sides.

  assert!(approx_equiv_v(projector.proj_int_mons_to_span_oshape_side_supp_basis_els([one, y, y*y], OShape(0), right_side),
                         [PolyBorrowingMons::new(~[1.,0.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,1.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,1.], vert_side_mons)], 1e-10));
  
  assert!(approx_equiv_v(projector.proj_int_mons_to_span_oshape_side_supp_basis_els([one, y, y*y], OShape(0), left_side),
                         [PolyBorrowingMons::new(~[1.,0.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,1.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,1.], vert_side_mons)], 1e-10));


  // Test projections with a perpendicular axis variable factor in the interior monomial, which should introduce a constant factor
  // into the projection which is 0 if projecting onto the lesser side along the perpendicular axis, or the side length along the
  // axis for the greater side. The fe side lengths are 2 horizontally and 3 vertically.
  
  assert!(approx_equiv_v(projector.proj_int_mons_to_span_oshape_side_supp_basis_els([x, x*y, x*y*y], OShape(0), right_side),
                         [PolyBorrowingMons::new(~[2.,0.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,2.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,2.], vert_side_mons)], 1e-10));
  
  assert!(approx_equiv_v(projector.proj_int_mons_to_span_oshape_side_supp_basis_els([x, x*y, x*y*y], OShape(0), left_side),
                         [PolyBorrowingMons::new(~[0.,0.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,0.], vert_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,0.], vert_side_mons)], 1e-10));
  
  assert!(approx_equiv_v(projector.proj_int_mons_to_span_oshape_side_supp_basis_els([y, y*x, y*x*x], OShape(0), top_side),
                         [PolyBorrowingMons::new(~[3.,0.,0.], horz_side_mons),
                          PolyBorrowingMons::new(~[0.,3.,0.], horz_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,3.], horz_side_mons)], 1e-10));
  
  assert!(approx_equiv_v(projector.proj_int_mons_to_span_oshape_side_supp_basis_els([y, y*x, y*x*x], OShape(0), bottom_side),
                         [PolyBorrowingMons::new(~[0.,0.,0.], horz_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,0.], horz_side_mons),
                          PolyBorrowingMons::new(~[0.,0.,0.], horz_side_mons)], 1e-10));
}


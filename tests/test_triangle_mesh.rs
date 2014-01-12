use common::{R, pow, Dim, Deg, DEFAULT_INTEGRATION_REL_ERR, DEFAULT_INTEGRATION_ABS_ERR};
use monomial::{Monomial, Mon2d};
use polynomial::{poly};
use vector_monomial::VectorMonomial;
use mesh::{Mesh, FENum, OShape, SideFace, NBSideNum, NBSideInclusions};
use triangle_mesh::*;
use triangle_mesh_builder::*;

use std::io::buffered::BufferedReader;
use std::io::mem::BufReader;
use std::num::{sqrt, abs};

static intg_tol: R = 1e-12;

// Test retrieval of side face endpoints.

// Return all side face endpoint pairs, with endpoints in traversal order.

#[test]
fn test_sideface_endpoint_pairs_vertex_order() {
  let (v0, v1, v2) = ((0.,0.), (1.,1.), (-1.,1.));
  assert_eq!(side_face_endpoint_pair(SideFace(0), v0,v1,v2, (2u8, 1u8, 2u8), VertexOrder), ((0.,0.),(0.5,0.5)));
  assert_eq!(side_face_endpoint_pair(SideFace(1), v0,v1,v2, (2u8, 1u8, 2u8), VertexOrder), ((0.5,0.5),(1.,1.)));
  assert_eq!(side_face_endpoint_pair(SideFace(2), v0,v1,v2, (2u8, 1u8, 2u8), VertexOrder), ((1.,1.),(-1.,1.)));
  assert_eq!(side_face_endpoint_pair(SideFace(3), v0,v1,v2, (2u8, 1u8, 2u8), VertexOrder), ((-1.,1.),(-0.5,0.5)));
  assert_eq!(side_face_endpoint_pair(SideFace(4), v0,v1,v2, (2u8, 1u8, 2u8), VertexOrder), ((-0.5,0.5), (0.,0.)));
}

// Same as above but order endpoints within endpoint pairs in a standardized order with lesser point first.
#[test]
fn test_sideface_endpoint_pairs_canon_order() {
  let (v0, v1, v2) = ((0.,0.), (1.,1.), (-1.,1.));
  assert_eq!(side_face_endpoint_pair(SideFace(0), v0,v1,v2, (2u8, 1u8, 2u8), LesserEndpointsFirst), ((0.,0.),(0.5,0.5)));
  assert_eq!(side_face_endpoint_pair(SideFace(1), v0,v1,v2, (2u8, 1u8, 2u8), LesserEndpointsFirst), ((0.5,0.5),(1.,1.)));
  assert_eq!(side_face_endpoint_pair(SideFace(2), v0,v1,v2, (2u8, 1u8, 2u8), LesserEndpointsFirst), ((-1.,1.),(1.,1.)));
  assert_eq!(side_face_endpoint_pair(SideFace(3), v0,v1,v2, (2u8, 1u8, 2u8), LesserEndpointsFirst), ((-1.,1.),(-0.5,0.5)));
  assert_eq!(side_face_endpoint_pair(SideFace(4), v0,v1,v2, (2u8, 1u8, 2u8), LesserEndpointsFirst), ((-0.5,0.5), (0.,0.)));
}

// Test RefTri calculation of normals.
#[test]
fn test_ref_tri_normals() {
  let ref_tri = RefTri::new((0.,0.), (1.,1.), (-1.,1.), 1., (1u8,2u8,3u8));
  assert_eq!(ref_tri.outward_normals_by_side_face,
    ~[(1./sqrt(2.),-1./sqrt(2.)),
      (0.,1.),
      (0.,1.),
      (-1./sqrt(2.),-1./sqrt(2.)),
      (-1./sqrt(2.),-1./sqrt(2.)),
      (-1./sqrt(2.),-1./sqrt(2.))]);
}

// Test integration between upper and lower lines converging to a point.

/*
 Integrate the one monomial between lines of slope -3 and 2, diverging from point (-1,-2) to x = 3.
 Horizontal span of generated triangle is 4.
 According to the slopes, the right face of the resulting triangle has length 3*4 + 2*4 = 20.
 Area of triangle = 1/2 * 20 * 4 = 40.
*/
#[test]
fn test_1_intg_btw_pt_and_vert_seg() {
  let one = |_: &[R]| { 1. };
  let max_var_deg = Deg(0);
  assert!(intg_poly_fn_btw_pt_and_vert_seg(one, max_var_deg, (-1.,-2.), -3., 2., 3.).approx_eq(&40.));
}

/*
 Monomial xy between lines diverging from point (2,1) of slope +/- 1, from x = 2 to x = 3.
 int_{x=2..3} int_{y=1+(x-2)..y=1-(x-2)} x y dy dx
  = int_{x=2..3} [1/2 x y^2]_|{y=1+(x-2)..1-(x-2)}
  = int_{x=2..3} x/2 ( (1+(x-2))^2 - (1-(x-2))^2 )  # (1+(x-2) + 1-(x-2)) (1+(x-2) - 1-(x-2)) = 2 (2(x-2)) = 4(x-2)
  = int_{x=2..3} x/2 4 (x-2)
  = int_{x=2..3} 2 x (x-2)
  = int_{x=2..3} 2x^2 - 4x
  = [2/3 x^3 - 2x^2]_|{x=2..3}
  = 2 + 2/3
*/
#[test]
fn test_xy_intg_btw_pt_and_vert_seg() {
  assert!(intg_poly_fn_btw_pt_and_vert_seg(|x| x[0]*x[1], Deg(1), (2.,1.), -1., 1., 3.)
            .approx_eq(&(2. + 2./3.)));
}

#[test]
fn test_xy_intg_btw_pt_and_vert_seg_flipped() {
  assert!(intg_poly_fn_btw_pt_and_vert_seg(|x| x[0]*x[1], Deg(1), (-2.,1.), -1., 1., -3.)
            .approx_eq(&(-(2. + 2./3.))));
}

// Do the same integrals this time with general function integration (adaptive quadrature).

#[test]
fn test_gen_fn_1_intg_btw_pt_and_vert_seg() {
  let one = |_: &[R]| { 1. };
  let max_var_deg = Deg(0);
  assert!(intg_fn_btw_pt_and_vert_seg(one, (-1.,-2.), -3., 2., 3., intg_tol, intg_tol)
            .approx_eq(&40.));
}

#[test]
fn test_gen_fn_xy_intg_btw_pt_and_vert_seg() {
  assert!(intg_fn_btw_pt_and_vert_seg(|x:&[R]| x[0]*x[1], (2.,1.), -1., 1., 3., intg_tol, intg_tol)
            .approx_eq(&(2. + 2./3.)));
}

#[test]
fn test_gen_fn_xy_intg_btw_pt_and_vert_seg_flipped() {
  assert!(intg_fn_btw_pt_and_vert_seg(|x:&[R]| x[0]*x[1], (-2.,1.), -1., 1., -3., intg_tol, intg_tol)
            .approx_eq(&(-(2. + 2./3.))));
}

/*
 * input mesh, before subdivision
 *    (0,3)
 *       | .
 *       |    .
 *  (0,0)|_______.(4,0)
 *  The base triangle's nodes are enumerated counterclockwise starting at (0,0).
 */

#[test]
fn test_right_tri_mesh_1_subdiv_basic_properties() {
  let msh_is = &mut str_rdr(msh_1);
  let mesh: TriMesh<Mon2d> = TriMeshBuilder::from_gmsh_msh_stream(msh_is, 1u /*subdiv iters*/,
                                                                  intg_tol, intg_tol,
                                                                  false); // don't load tags
  assert_eq!(mesh.num_fes(), 4);
  assert_eq!(mesh.num_nb_sides(), 3);
  assert_eq!(mesh.num_oriented_element_shapes(), 2);
  assert_eq!(mesh.max_fe_diameter(), 2.5);
}

#[test]
fn test_right_tri_mesh_1_subdiv_subtris() {
  let msh_is = &mut str_rdr(msh_1);
  let mesh: TriMesh<Mon2d> = TriMeshBuilder::from_gmsh_msh_stream(msh_is, 1u /*subdiv iters*/,
                                                                  intg_tol, intg_tol,
                                                                  false); // don't load tags
  // 3 primary subtriangles
  assert_eq!(mesh.oriented_shape_for_fe(FENum(0)), OShape(0));
  assert_eq!(mesh.fes[0].oshape, OShape(0));
  assert_eq!(mesh.fes[0].v0, (0.,0.));
  assert_eq!(mesh.oriented_shape_for_fe(FENum(1)), OShape(0));
  assert_eq!(mesh.fes[1].oshape, OShape(0));
  assert_eq!(mesh.fes[1].v0, (2.,0.));
  assert_eq!(mesh.oriented_shape_for_fe(FENum(2)), OShape(0));
  assert_eq!(mesh.fes[2].oshape, OShape(0));
  assert_eq!(mesh.fes[2].v0, (0.,3./2.));

  // secondary subtriangle
  assert_eq!(mesh.oriented_shape_for_fe(FENum(3)), OShape(1));
  assert_eq!(mesh.fes[3].oshape, OShape(1));
  assert_eq!(mesh.fes[3].v0, (2.,0.));
}

#[test]
fn test_right_tri_mesh_1_subdiv_primary_oshape() {
  let msh_is = &mut str_rdr(msh_1);
  let mesh: TriMesh<Mon2d> = TriMeshBuilder::from_gmsh_msh_stream(msh_is, 1u /*subdiv iters*/,
                                                                  intg_tol, intg_tol,
                                                                  false); // don't load tags
  assert_eq!(mesh.num_side_faces_for_oshape(OShape(0)), 3);
  assert_eq!(mesh.oshapes[0].v01, (2.,0.));
  assert_eq!(mesh.oshapes[0].v02, (0.,3./2.));
  assert_eq!(mesh.oshapes[0].outward_normals_by_side_face, ~[(0.,-1.), (3./5.,4./5.), (-1.,-0.)]);
  assert_eq!(mesh.shape_diameter_inv(OShape(0)), 1./2.5);
  assert_eq!(mesh.dependent_dim_for_oshape_side(OShape(0), SideFace(0)), Dim(1));
  assert_eq!(mesh.dependent_dim_for_oshape_side(OShape(0), SideFace(1)), Dim(1));
  assert_eq!(mesh.dependent_dim_for_oshape_side(OShape(0), SideFace(2)), Dim(0));
}

#[test]
fn test_right_tri_mesh_1_subdiv_secondary_oshape() {
  let msh_is = &mut str_rdr(msh_1);
  let mesh: TriMesh<Mon2d> = TriMeshBuilder::from_gmsh_msh_stream(msh_is, 1u /*subdiv iters*/,
                                                                  intg_tol, intg_tol,
                                                                  false); // don't load tags
  assert_eq!(mesh.num_side_faces_for_oshape(OShape(1)), 3);
  assert_eq!(mesh.oshapes[1].v01, (0.,3./2.));
  assert_eq!(mesh.oshapes[1].v02, (-2.,3./2.));
  assert_eq!(mesh.oshapes[1].outward_normals_by_side_face, ~[(1.,-0.), (0.,1.), (-3./5.,-4./5.)]);
  assert_eq!(mesh.shape_diameter_inv(OShape(1)), 1./2.5);
  assert_eq!(mesh.dependent_dim_for_oshape_side(OShape(1), SideFace(0)), Dim(0));
  assert_eq!(mesh.dependent_dim_for_oshape_side(OShape(1), SideFace(1)), Dim(1));
  assert_eq!(mesh.dependent_dim_for_oshape_side(OShape(1), SideFace(2)), Dim(1));
}

#[test]
fn test_right_tri_mesh_1_subdiv_boundary_sides() {
  let msh_is = &mut str_rdr(msh_1);
  let mesh: TriMesh<Mon2d> = TriMeshBuilder::from_gmsh_msh_stream(msh_is, 1u /*subdiv iters*/,
                                                                  intg_tol, intg_tol,
                                                                  false); // don't load tags
  assert_eq!(mesh.num_boundary_sides(), 6);
  assert!( mesh.is_boundary_side(FENum(0), SideFace(0)));
  assert!(!mesh.is_boundary_side(FENum(0), SideFace(1)));
  assert!( mesh.is_boundary_side(FENum(0), SideFace(2)));
  assert!( mesh.is_boundary_side(FENum(1), SideFace(0)));
  assert!( mesh.is_boundary_side(FENum(1), SideFace(1)));
  assert!(!mesh.is_boundary_side(FENum(1), SideFace(2)));
  assert!(!mesh.is_boundary_side(FENum(2), SideFace(0)));
  assert!( mesh.is_boundary_side(FENum(2), SideFace(1)));
  assert!( mesh.is_boundary_side(FENum(2), SideFace(2)));
}

#[test]
fn test_right_tri_mesh_1_subdiv_nb_side_incls() {
  let msh_is = &mut str_rdr(msh_1);
  let mesh: TriMesh<Mon2d> = TriMeshBuilder::from_gmsh_msh_stream(msh_is, 1u /*subdiv iters*/,
                                                                  intg_tol, intg_tol,
                                                                  false); // don't load tags
  assert_eq!(mesh.nbsideincls_by_nbsidenum.len(), 3);
  let nbsn_fe0sf1 = mesh.nb_side_num_for_fe_side(FENum(0), SideFace(1));
  assert_eq!(mesh.nb_side_num_for_fe_side(FENum(3), SideFace(2)), nbsn_fe0sf1);
  assert_eq!(mesh.fe_inclusions_of_nb_side(nbsn_fe0sf1),
             NBSideInclusions { nb_side_num: nbsn_fe0sf1,
                                fe1: FENum(0), side_face_in_fe1: SideFace(1),
                                fe2: FENum(3), side_face_in_fe2: SideFace(2)});
  
  assert_eq!(mesh.nbsideincls_by_nbsidenum.len(), 3);
  let nbsn_fe1sf2 = mesh.nb_side_num_for_fe_side(FENum(1), SideFace(2));
  assert_eq!(mesh.nb_side_num_for_fe_side(FENum(3), SideFace(0)), nbsn_fe1sf2);
  assert_eq!(mesh.fe_inclusions_of_nb_side(nbsn_fe1sf2),
             NBSideInclusions { nb_side_num: nbsn_fe1sf2,
                                fe1: FENum(1), side_face_in_fe1: SideFace(2),
                                fe2: FENum(3), side_face_in_fe2: SideFace(0)});
  
  assert_eq!(mesh.nbsideincls_by_nbsidenum.len(), 3);
  let nbsn_fe2sf0 = mesh.nb_side_num_for_fe_side(FENum(2), SideFace(0));
  assert_eq!(mesh.nb_side_num_for_fe_side(FENum(3), SideFace(1)), nbsn_fe2sf0);
  assert_eq!(mesh.fe_inclusions_of_nb_side(nbsn_fe2sf0),
             NBSideInclusions { nb_side_num: nbsn_fe2sf0,
                                fe1: FENum(2), side_face_in_fe1: SideFace(0),
                                fe2: FENum(3), side_face_in_fe2: SideFace(1)});
}

/*
# Test integrals.

# primary reference element
#     (0,1.5)
#       | .
#       |    .
#  (0,0)|_______.(2,0)
#

# integrals of face-relative monomials on side faces

# integrals of constant 1 monomials should give side lengths
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(1), tmsh),
               2., atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(2), tmsh),
               hypot(2.,1.5), rtol=1e-15, atol=1e-12)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(3), tmsh),
               1.5, rtol=1e-15, atol=1e-12)

@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)

# These compute side lengths like the above, but are expressed as a function times a local monomial.
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(1), tmsh),
               2., atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(2), tmsh),
               hypot(2.,1.5), atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(3), tmsh),
               1.5, atol=1e-15, rtol=1e-15)

# integral of xy on side 1
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y, fenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)
xy_on_fe1_side1 = local_mon_on_fe_side_as_global_fn(deg(1), deg(1), fenum(1), fefacenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side1, onemon, fenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)
x_on_fe1_side1 = local_mon_on_fe_side_as_global_fn(deg(1), deg(0), fenum(1), fefacenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(x_on_fe1_side1, y, fenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)

# Integrate polynomial xy on side face 2 of fe/oshape 1.
# We can traverse side 2 of fe 1's reference triangle in side relative coordinates with
# p(t) = (2,0)+t((0,1.5)-(2,0))-(1,0.75), 0<=t<=1. Then
# xy(p(t)) = xy((1 - 2t, 1.5t - 0.75)) = (1-2t)(1.5t-0.75)
# so the integral should be
# int_{0,..1} (1-2t)(1.5t-0.75) |(-2,1.5)| dt = 2.5[-.75t + 1.5t^2 - t^3]|{t=0,1} = -0.625
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y, fenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)
x_on_fe1_side2 = local_mon_on_fe_side_as_global_fn(deg(1), deg(0), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(x_on_fe1_side2, y, fenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)
xy_on_fe1_side2 = local_mon_on_fe_side_as_global_fn(deg(1), deg(1), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side2, onemon, fenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)

# Integrate the same monomial over the same side but now also vs. outward normal.
one_comp1_vmon = VectorMonomial(onemon, dim(1))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x*y, one_comp1_vmon, oshapenum(1), fefacenum(2), tmsh),
               -0.625 * 3/5, atol=1e-15, rtol=1e-15) # 3/5 = component 1 of outward normal
one_comp2_vmon = VectorMonomial(onemon, dim(2))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x*y, one_comp2_vmon, oshapenum(1), fefacenum(2), tmsh),
               -0.625 * 4/5, atol=1e-15, rtol=1e-15) # 4/5 = component 2 of outward normalj
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(onemon, one_comp1_vmon, oshapenum(1), fefacenum(1), tmsh),
               0, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(onemon, one_comp2_vmon, oshapenum(1), fefacenum(3), tmsh),
               0, atol=1e-15, rtol=1e-15)

# Integrate (0,x) vector monomial vs outward normal along first side, which should equal int_0^2 -x dx = -2.
x_comp2_vmon = VectorMonomial(x, dim(2))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(onemon, x_comp2_vmon, oshapenum(1), fefacenum(1), tmsh),
               -2, atol=1e-15, rtol=1e-15)

# Like the above but with a scalar monomial instead of vector dotted against the normal.
@test isapprox(Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x, onemon, oshapenum(1), fefacenum(1), tmsh),
               2, atol=1e-15, rtol=1e-15)

# Integrate x side-local monomial vs (0,x) vector monomial (fe-relative) vs outward normal along the first side.
# As an fe-local monomial, the side local monomial x is x-1. Thus the integral should be
# int_0^2 -(x-1)x dx = -1/3 2^3 + 2^2/2 = -2/3
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x, x_comp2_vmon, oshapenum(1), fefacenum(1), tmsh),
               -2/3, atol=1e-15, rtol=1e-15)
# Like the above but with a scalar monomial instead of vector dotted against the normal.
@test isapprox(Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x, x, oshapenum(1), fefacenum(1), tmsh),
               2/3, atol=1e-15, rtol=1e-15)

# Integrate y side-local monomial vs (y,0) fe-relative vector monomial vs outward normal along the third side.
# As an fe-local monomial, the side local monomial y is y-3/4. Thus the integral should be
# int_0^{3/2} -(y-3/4)y dy = -1/2 (3/2)^2 + 3/8 (3/2)^2 = -9/32
y_comp1_vmon = VectorMonomial(y, dim(1))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y, y_comp1_vmon, oshapenum(1), fefacenum(3), tmsh),
               -9/32, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y, y, oshapenum(1), fefacenum(3), tmsh),
               9/32, atol=1e-15, rtol=1e-15)


# interior integrals

# Integral of constant one monomial on a triangle interior should give its area.
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(0), tmsh),
               0.5*2*1.5, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(2), fefacenum(0), tmsh),
               0.5*2*1.5, atol=1e-15, rtol=1e-15)

# integral xy^2 over the interior of reference triangle 1
# = int_{x=0..2} int_{y=0..1.5-0.75x} xy^2 = int_{x=0..2} 1/3 x(1.5-0.75x)^3
# = 1/3 (1.6875 x^2 + -1.6875 x^3 + 0.6328125 x^4 + -0.084375 x^5)|_{x=0,2}
# = 0.225
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y^2, oshapenum(1), fefacenum(0), tmsh),
               0.225, atol=1e-15, rtol=1e-15)


# Integrate function (x,y)->(x,y)-o(fe) against local monomial y on various finite element interiors, which
# should equal the integral of xy^2 local monomial over the oriented shape's interior.
xy_on_fe1 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1),fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(0), tmsh),
               0.225, rtol=1e-13, atol=1e-13)
xy_on_fe2 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1),fenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe2, y, fenum(2), fefacenum(0), tmsh),
               0.225, rtol=1e-13, atol=1e-13)
xy_on_fe3 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1),fenum(3), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe3, y, fenum(3), fefacenum(0), tmsh),
               0.225, rtol=1e-13, atol=1e-13)

# input mesh, before subdivision
#     (2,8).
#         . .
#        .   .
#  (0,0). . . .(4,0)
#
isosc_tri_mesh =
    """
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$Nodes
    3
    1 0 0 0
    2 4 0 0
    3 2 8 0
    \$EndNodes
    \$Elements
    7
    1 15 2 0 1 1
    2 15 2 0 2 2
    3 15 2 0 3 3
    4 1 2 0 1 1 2
    5 1 2 0 2 2 3
    6 1 2 0 3 3 1
    7 2 2 0 5 1 2 3
    \$EndElements
    """

# Check that integrals that must be done in two pieces are done properly.
tmsh = TriMesh(IOString(isosc_tri_mesh), 2)
@test tmsh.fes[1] == ElTri(oshapenum(1),Vec(0.,0.))

# primary reference element
#   (1/2,2).
#         . .
#        .   .
#  (0,0). . . .(1,0)
#

@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(2), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)

# Integrate the local monomial xy^2 over the interior of a reference triangle.
# int_0^2 int_{y/4}^{-y/4 + 1}  xy^2 dx dy
#  = int_0^2 1/2 y^2 ((-y/4 + 1)^2 - (y/4)^2) dy
#  = int_0^2 1/2 y^2 - 1/4 y^3 dy
#  = (1/6 y^3 - 1/16 y^4)|_0^2
#  = 8/6 - 1 = 1/3
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x * y^2, oshapenum(1), fefacenum(0), tmsh),
               1/3, rtol=1e-11, atol=1e-11)
# Like the above, but expressed as a product of a global function, chosen to match local monomial xy, and monomial y.
xy_on_fe1 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1), fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(0), tmsh),
               1/3, rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)

# Integrate xy^2 monomial along side face 1, on which y = 0.
xy_on_fe1 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1), fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(1), tmsh),
               0., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y^2, fenum(1), fefacenum(1), tmsh),
               0., rtol=1e-11, atol=1e-11)

# Integrate polynomial xy on side face 2 of fe/oshape 1.
# We can traverse side 2 of fe 1's reference triangle in side relative coordinates with
# p(t) = (1,0) + t((1/2,2)-(1,0)) - (0.75,1)
#      = (1/4,-1) + t(-1/2,2),   (0<=t<=1). Then
# xy(p(t)) = xy((1/4 - 1/2 t, -1 + 2t)) = -1/4 + t - t^2
# so the integral should be
# int_0^1 (-1/4 + t - t^2) |(-1/2,2)| dt = sqrt(17)/2 [-t/4 + 1/2 t^2 -1/3 t^3]|_{t=0,1} = -sqrt(17)/24
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y, fenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)
x_on_fe1_side2 = local_mon_on_fe_side_as_global_fn(deg(1), deg(0), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(x_on_fe1_side2, y, fenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)
xy_on_fe1_side2 = local_mon_on_fe_side_as_global_fn(deg(1), deg(1), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side2, onemon, fenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)


# inversion of the above, original undivided triangle
#  (-2,8). . . .(2,8)
#         .   .
#          . .
#           .(0,0)
inv_isosc_tri_mesh =
    """
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$Nodes
    3
    1 0 0 0
    2 2 8 0
    3 -2 8 0
    \$EndNodes
    \$Elements
    7
    1 15 2 0 1 1
    2 15 2 0 2 2
    3 15 2 0 3 3
    4 1 2 0 1 1 2
    5 1 2 0 2 2 3
    6 1 2 0 3 3 1
    7 2 2 0 5 1 2 3
    \$EndElements
    """

# primary reference element
# (-1/2,2). . . .(1/2,2)
#          .   .
#           . .
#            .(0,0)

# Check that integrals that must be done in two pieces are done properly.
tmsh = TriMesh(IOString(inv_isosc_tri_mesh), 2)
@test tmsh.fes[1] == ElTri(oshapenum(1),Vec(0.,0.))

@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(2), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)

# Integrate the local monomial xy^2 over the interior of a reference triangle.
# int_0^2 int_{-y/4}^{y/4}  xy^2 dx dy
#  = int_0^2 1/2 y^2 ((y/4)^2 - (-y/4)^2) dy
#  = 0
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x * y^2, oshapenum(1), fefacenum(0), tmsh),
               0., rtol=1e-11, atol=1e-11)

# Like the above, but expressed as a product of a global function, chosen to match local monomial xy, and monomial y.
xy_on_fe1 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1), fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(0), tmsh),
               0., rtol=1e-11, rtol=1e-11)

@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)

# Integrate the local monomial x^2y^2 over the interior of the primary reference triangle.
# int_0^2 int_{-y/4}^{y/4}  x^2y^2 dx dy
#  = int_0^2 1/3 y^2 ((y/4)^3 - (-y/4)^3) dy
#  = int_0^2 1/3 y^2 2/64 y^3 dy
#  = int_0^2 1/96 y^5 dy
#  = 1/96 1/6 y^6|_0^2
#  = 1/9
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x^2 * y^2, oshapenum(1), fefacenum(0), tmsh),
               1/9, rtol=1e-11, atol=1e-11)

# Same as above, but expressed as a product of a global function and monomial.

@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x^2 * y^2, fenum(1), fefacenum(0), tmsh),
               1/9, rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(local_mon_on_fe_int_as_global_fn(deg(2),deg(2), fenum(1), tmsh), onemon, fenum(1), fefacenum(0), tmsh),
               1/9, rtol=1e-11, atol=1e-11)

# areas
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(2), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(3), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(4), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)

# Integrate polynomial x^3y on side face 2 of fe/oshape 1.
# The integral is int_{-1/2}^{1/2} 2 x^3 dx = 0.
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x^3*y, oshapenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x^3*y, fenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)
xcubed_on_fe1_side2 = local_mon_on_fe_side_as_global_fn(deg(3), deg(0), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xcubed_on_fe1_side2, y, fenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)
xy_on_fe1_side2 = local_mon_on_fe_side_as_global_fn(deg(1), deg(1), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side2, x^2, fenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)


# secondary reference element
#  (-1/2,2).
#         . .
#        .   .
# (-1,0). . . .(0,0)

# Integrate the local monomial xy^2 over the interior of the secondary reference triangle.
# int_0^2 int_{y/4-1}^{-y/4}  xy^2 dx dy
#  = int_0^2 1/2 y^2 ((-y/4)^2 - (y/4-1)^2) dy
#  = int_0^2 1/4 y^3 - 1/2 y^2 dy
#  = (1/16 y^4 - 1/6 y^3)|_0^2
#  = 1 - 8/6 = -1/3
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x * y^2, oshapenum(2), fefacenum(0), tmsh),
               -1/3, atol=1e-11, rtol=1e-11)

# Like the above, but expressed as a product of a global function, chosen to match local monomial xy, and monomial y.
xy_on_fe4 = local_mon_on_fe_int_as_global_fn(deg(1),deg(1), fenum(4), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe4, y, fenum(4), fefacenum(0), tmsh),
               -1/3, atol=1e-11, rtol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y^2, fenum(4), fefacenum(0), tmsh),
               -1/3, atol=1e-11, rtol=1e-11)


# Test construction from rectangle mesh.
rmesh2x2 = RectMesh([1.,1.], [3.,3.], [RMesh.mesh_coord(2), RMesh.mesh_coord(2)])
tmesh = from_rect_mesh(rmesh2x2, 1)
@test Mesh.num_fes(tmesh) == 32
@test Mesh.num_boundary_sides(tmesh) == 16
*/

#[inline]
pub fn midpt(p: Point, q: Point) -> Point { (0.5*(p.n0() + q.n0()), 0.5*(p.n1() + q.n1())) }

#[inline]
fn vsum(p: Vec, q: Vec) -> Vec { (p.n0() + q.n0(), p.n1() + q.n1()) }


// Local monomials extended as global functions.

fn local_mon_on_fe_int_as_global_fn <Mon: Monomial>
   ( exp1: Deg,
     exp2: Deg, 
     fe: FENum,
     mesh: &TriMesh<Mon> )
   -> proc(&[R]) -> R
{
  let o = mesh.el_tri(fe).v0;
  proc(x: &[R]) { pow(x[0] - o.n0(), *exp1 as uint) * pow(x[1] - o.n1(), *exp2 as uint) }
}


fn local_mon_on_fe_side_as_global_fn <Mon: Monomial>
   ( exp1: Deg,
     exp2: Deg, 
     fe: FENum,
     sf: SideFace,
     mesh: &TriMesh<Mon> )
   -> proc(&[R]) -> R
{
  let ref_tri = mesh.ref_tri(fe);
  let v0 = mesh.el_tri(fe).v0;
  let (a,b) = side_face_endpoint_pair(sf, v0, vsum(v0, ref_tri.v01), vsum(v0, ref_tri.v02),
                                      ref_tri.nums_side_faces_between_vertexes, VertexOrder);
  let o = midpt(a, b); // The local origin of each side face is the midpoint.
  proc(x: &[R]) { pow(x[0] - o.n0(), *exp1 as uint) * pow(x[1] - o.n1(), *exp2 as uint) }
}

fn str_rdr<'a>(s: &'a str) -> BufferedReader<BufReader<'a>> {
  BufferedReader::new(BufReader::new(s.as_bytes()))
}

static msh_1: &'static str = 
"$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0 0 0
2 4 0 0
3 0 3 0
$EndNodes
$Elements
7
1 15 2 0 1 1
2 15 2 0 2 2
3 15 2 0 3 3
4 1 2 0 1 1 2
5 1 2 0 2 2 3
6 1 2 0 3 3 1
7 2 2 0 5 1 2 3
$EndElements
";


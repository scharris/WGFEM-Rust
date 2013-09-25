use std::num::{sqrt, abs};
use common::*;
use monomial::{Mon1d, Mon2d, Mon3d, Mon4d};
use polynomial::{poly};
use vector_monomial::VectorMonomial;
use mesh::*;
use rectangle_mesh::*;

#[test]
fn test_3x4_constr() -> () {
  let mesh_min_coords = ~[1f64, 2.];
  let mesh_max_coords = ~[2f64, 3.];
  let mesh_ldims = ~[MeshCoord(3), MeshCoord(4)];
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new_with_err_tols(mesh_min_coords.clone(), mesh_max_coords.clone(), mesh_ldims.clone(), 1e-4, 1e-5);

  assert_eq!(rmesh3x4.space_dims, 2);
  assert_eq!(&rmesh3x4.min_bounds, &mesh_min_coords);
  assert_eq!(&rmesh3x4.max_bounds, &mesh_max_coords);
  assert_eq!(&rmesh3x4.mesh_ldims, &mesh_ldims);
  assert_eq!(&rmesh3x4.fe_side_lens, &~[1./3., 1./4.]);
  assert_approx(rmesh3x4.rect_diameter, sqrt(pow(1./3.,2) + pow(1./4.,2)));
  assert_approx(rmesh3x4.shape_diameter_inv(OShape(0)), 1./sqrt(pow(1./3.,2) + pow(1./4.,2)));
  assert_eq!(&rmesh3x4.cumprods_mesh_ldims, &~[3, 3*4]);

  assert_eq!(&rmesh3x4.cumprods_nb_side_mesh_ldims_by_perp_axis, &~[
             ~[3-1, (3-1)*4],
             ~[3, 3*(4-1)]]);

  assert_eq!(&rmesh3x4.first_nb_side_nums_by_perp_axis, &~[NBSideNum(0u), NBSideNum((3-1)*4),]);

  assert_eq!(rmesh3x4.num_fes(), 3*4);
  assert_eq!(rmesh3x4.num_nb_sides(), (3-1)*4 + 3*(4-1));
  assert_eq!(rmesh3x4.num_side_faces_per_fe, 4);
  assert_eq!(rmesh3x4.num_oriented_element_shapes(), 1);

  assert_eq!(rmesh3x4.integration_rel_err, 1e-4);
  assert_eq!(rmesh3x4.integration_abs_err, 1e-5);
}

#[test]
fn test_3x4x5_constr() -> () {
  let mesh_min_coords = ~[1f64, 2., 3.];
  let mesh_max_coords = ~[2f64, 3., 4.];
  let mesh_ldims = ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)];
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(mesh_min_coords.clone(), mesh_max_coords.clone(), mesh_ldims.clone());

  assert_eq!(rmesh3x4x5.space_dims, 3);
  assert_eq!(&rmesh3x4x5.min_bounds, &mesh_min_coords);
  assert_eq!(&rmesh3x4x5.max_bounds, &mesh_max_coords);
  assert_eq!(&rmesh3x4x5.mesh_ldims, &mesh_ldims);
  assert_eq!(&rmesh3x4x5.fe_side_lens, &~[1./3., 1./4., 1./5.]);
  assert_approx(rmesh3x4x5.rect_diameter, sqrt(pow(1./3.,2) + pow(1./4.,2) + pow(1./5.,2)));
  assert_approx(rmesh3x4x5.shape_diameter_inv(OShape(0)), 1./sqrt(pow(1./3.,2) + pow(1./4.,2) + pow(1./5.,2)));
  assert_eq!(&rmesh3x4x5.cumprods_mesh_ldims, &~[3, 3*4, 3*4*5]);

  assert_eq!(&rmesh3x4x5.cumprods_nb_side_mesh_ldims_by_perp_axis, &~[
             ~[3-1, (3-1)*4, (3-1)*4*5],
             ~[3, 3*(4-1), 3*(4-1)*5],
             ~[3, 3*4, 3*4*(5-1)]]);

  assert_eq!(&rmesh3x4x5.first_nb_side_nums_by_perp_axis, &~[NBSideNum(0u), NBSideNum((3-1)*4*5), NBSideNum((3-1)*4*5 + 3*(4-1)*5)]);

  assert_eq!(rmesh3x4x5.num_fes(), 3*4*5);
  assert_eq!(rmesh3x4x5.num_nb_sides(), (3-1)*4*5 + 3*(4-1)*5 + 3u * 4 * (5-1));
  assert_eq!(rmesh3x4x5.num_side_faces_per_fe, 6);
  assert_eq!(rmesh3x4x5.num_oriented_element_shapes(), 1);

  assert_eq!(rmesh3x4x5.integration_rel_err, DEFAULT_INTEGRATION_REL_ERR);
  assert_eq!(rmesh3x4x5.integration_abs_err, DEFAULT_INTEGRATION_ABS_ERR);
}

#[test]
fn test_3x4x5x6_constr() -> () {
  let mesh_min_coords = ~[1f64, 2., 3., 4.];
  let mesh_max_coords = ~[2f64, 3., 4., 5.];
  let mesh_ldims = ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)];
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(mesh_min_coords.clone(), mesh_max_coords.clone(), mesh_ldims.clone());

  assert_eq!(rmesh3x4x5x6.space_dims, 4);
  assert_eq!(&rmesh3x4x5x6.min_bounds, &mesh_min_coords);
  assert_eq!(&rmesh3x4x5x6.max_bounds, &mesh_max_coords);
  assert_eq!(&rmesh3x4x5x6.mesh_ldims, &mesh_ldims);
  assert_eq!(&rmesh3x4x5x6.fe_side_lens, &~[1./3., 1./4., 1./5., 1./6.]);
  assert_approx(rmesh3x4x5x6.rect_diameter, sqrt(pow(1./3.,2) + pow(1./4.,2) + pow(1./5.,2) + pow(1./6.,2)));
  assert_approx(rmesh3x4x5x6.shape_diameter_inv(OShape(0)), 1./sqrt(pow(1./3.,2) + pow(1./4.,2) + pow(1./5.,2) + pow(1./6.,2)));
  assert_eq!(&rmesh3x4x5x6.cumprods_mesh_ldims, &~[3, 3*4, 3*4*5, 3*4*5*6]);

  assert_eq!(&rmesh3x4x5x6.cumprods_nb_side_mesh_ldims_by_perp_axis, &~[
             ~[3-1, (3-1)*4, (3-1)*4*5, (3-1)*4*5*6],
             ~[3, 3*(4-1), 3*(4-1)*5, 3*(4-1)*5*6],
             ~[3, 3*4, 3*4*(5-1), 3*4*(5-1)*6],
             ~[3, 3*4, 3*4*5, 3*4*5*(6-1)]]);

  assert_eq!(&rmesh3x4x5x6.first_nb_side_nums_by_perp_axis,
             &~[NBSideNum(0u),
                NBSideNum((3-1)*4*5*6),
                NBSideNum((3-1)*4*5*6 + 3*(4-1)*5*6),
                NBSideNum((3-1)*4*5*6 + 3*(4-1)*5*6 + 3*4*(5-1)*6)]);

  assert_eq!(rmesh3x4x5x6.num_fes(), 3*4*5*6);
  assert_eq!(rmesh3x4x5x6.num_nb_sides(), (3-1)*4*5*6 + 3*(4-1)*5*6 + 3*4*(5-1)*6 + 3*4*5*(6-1));
  assert_eq!(rmesh3x4x5x6.num_side_faces_per_fe, 8);
  assert_eq!(rmesh3x4x5x6.num_oriented_element_shapes(), 1);

  assert_eq!(rmesh3x4x5x6.integration_rel_err, DEFAULT_INTEGRATION_REL_ERR);
  assert_eq!(rmesh3x4x5x6.integration_abs_err, DEFAULT_INTEGRATION_ABS_ERR);
}

#[test]
fn test_3x4_mesh_coords() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(0)), MeshCoord(0));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(0)), MeshCoord(0));
  
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(1)), MeshCoord(1));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(1)), MeshCoord(0));
  
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(2)), MeshCoord(2));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(2)), MeshCoord(0));

  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(3)), MeshCoord(0));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(3)), MeshCoord(1));

  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(5)), MeshCoord(2));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(5)), MeshCoord(1));

  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(6)), MeshCoord(0));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(6)), MeshCoord(2));

  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(9)), MeshCoord(0));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(9)), MeshCoord(3));

  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(0), FENum(11)), MeshCoord(2));
  assert_eq!(rmesh3x4.fe_mesh_coord(Dim(1), FENum(11)), MeshCoord(3));
}

#[test]
fn test_3x4x5_mesh_coords() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(0)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(0)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(0)), MeshCoord(0));
  
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(1)), MeshCoord(1));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(1)), MeshCoord(0));
  
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(2)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(2)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(2)), MeshCoord(0));

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(3)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(3)), MeshCoord(1));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(3)), MeshCoord(0));

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(5)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(5)), MeshCoord(1));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(5)), MeshCoord(0));

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(6)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(6)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(6)), MeshCoord(0));

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(9)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(9)), MeshCoord(3));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(9)), MeshCoord(0));

  // last element of first stack
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(11)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(11)), MeshCoord(3));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(11)), MeshCoord(0));

  // first element of second stack
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(12)), MeshCoord(1));

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(14)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(14)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(14)), MeshCoord(1));

  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(15)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(15)), MeshCoord(1));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(15)), MeshCoord(1));

  // last element of second stack
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(13+12-2)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(13+12-2)), MeshCoord(3));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(13+12-2)), MeshCoord(1));


  // first element of third stack
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(13+12-1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(13+12-1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(13+12-1)), MeshCoord(2));

  // first element of fifth and last stack
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(4*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(4*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(4*12)), MeshCoord(4));

  // last element of fifth and last stack
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(5*12-1)), MeshCoord(2));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(1), FENum(5*12-1)), MeshCoord(3));
  assert_eq!(rmesh3x4x5.fe_mesh_coord(Dim(2), FENum(5*12-1)), MeshCoord(4));
}

#[test]
fn test_3x4x5x6_mesh_coords() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(0)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(0)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(0)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(0)), MeshCoord(0));
  
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(1)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(1)), MeshCoord(0));
  
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(2)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(2)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(2)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(2)), MeshCoord(0));

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(3)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(3)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(3)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(3)), MeshCoord(0));

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(5)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(5)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(5)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(5)), MeshCoord(0));

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(6)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(6)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(6)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(6)), MeshCoord(0));

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(9)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(9)), MeshCoord(3));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(9)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(9)), MeshCoord(0));

  // last element of first 2-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(11)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(11)), MeshCoord(3));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(11)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(11)), MeshCoord(0));

  // first element of second 2-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(12)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(12)), MeshCoord(0));

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(14)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(14)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(14)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(14)), MeshCoord(0));

  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(15)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(15)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(15)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(15)), MeshCoord(0));

  // last element of second 2-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(13+12-2)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(13+12-2)), MeshCoord(3));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(13+12-2)), MeshCoord(1));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(13+12-2)), MeshCoord(0));


  // first element of third 2-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(13+12-1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(13+12-1)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(13+12-1)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(13+12-1)), MeshCoord(0));

  // first element of fifth and last 2-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(4*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(4*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(4*12)), MeshCoord(4));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(4*12)), MeshCoord(0));

  // last element of fifth and last 2-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(5*12-1)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(5*12-1)), MeshCoord(3));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(5*12-1)), MeshCoord(4));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(5*12-1)), MeshCoord(0));

  // first element of second 3-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(5*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(5*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(5*12)), MeshCoord(0));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(5*12)), MeshCoord(1));
  
  // last element of sixth and last 3-stack
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(3*4*5*6-1)), MeshCoord(2));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(1), FENum(3*4*5*6-1)), MeshCoord(3));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(2), FENum(3*4*5*6-1)), MeshCoord(4));
  assert_eq!(rmesh3x4x5x6.fe_mesh_coord(Dim(3), FENum(3*4*5*6-1)), MeshCoord(5));
}

#[test]
#[should_fail]
fn test_3x4_bad_mesh_coords() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  rmesh3x4.fe_mesh_coord(Dim(0), FENum(3*4));
}

#[test]
#[should_fail]
fn test_3x4x5_bad_mesh_coords() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(3*4*5));
}

#[test]
#[should_fail]
fn test_3x4x5x6_bad_mesh_coords() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  rmesh3x4x5x6.fe_mesh_coord(Dim(0), FENum(3*4*5*6));
}

fn test_3x4_boundary_side_determ() -> () {
  let rmesh3x4: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));

  assert!( rmesh3x4.is_boundary_side(FENum(0), left_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(0), right_side));
  assert!( rmesh3x4.is_boundary_side(FENum(0), bottom_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(0), top_side));

  assert!(!rmesh3x4.is_boundary_side(FENum(2), left_side));
  assert!( rmesh3x4.is_boundary_side(FENum(2), right_side));
  assert!( rmesh3x4.is_boundary_side(FENum(2), bottom_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(2), top_side));

  assert!( rmesh3x4.is_boundary_side(FENum(3), left_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(3), right_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(3), bottom_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(3), top_side));

  assert!(!rmesh3x4.is_boundary_side(FENum(11), left_side));
  assert!( rmesh3x4.is_boundary_side(FENum(11), right_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(11), bottom_side));
  assert!( rmesh3x4.is_boundary_side(FENum(11), top_side));

  assert!( rmesh3x4.is_boundary_side(FENum(12), left_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(12), right_side));
  assert!( rmesh3x4.is_boundary_side(FENum(12), bottom_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(12), top_side));

  assert!( rmesh3x4.is_boundary_side(FENum(4*12), left_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(4*12), right_side));
  assert!( rmesh3x4.is_boundary_side(FENum(4*12), bottom_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(4*12), top_side));

  assert!(!rmesh3x4.is_boundary_side(FENum(5*12-1), left_side));
  assert!( rmesh3x4.is_boundary_side(FENum(5*12-1), right_side));
  assert!(!rmesh3x4.is_boundary_side(FENum(5*12-1), bottom_side));
  assert!( rmesh3x4.is_boundary_side(FENum(5*12-1), top_side));
}


fn test_3x4x5_boundary_side_determ() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  assert!( rmesh3x4x5.is_boundary_side(FENum(0), left_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(0), right_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(0), bottom_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(0), top_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(0), back_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(0), front_side));

  assert!(!rmesh3x4x5.is_boundary_side(FENum(2), left_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(2), right_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(2), bottom_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(2), top_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(2), back_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(2), front_side));

  assert!( rmesh3x4x5.is_boundary_side(FENum(3), left_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), right_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), bottom_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), top_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(3), back_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), front_side));

  assert!(!rmesh3x4x5.is_boundary_side(FENum(11), left_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(11), right_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(11), bottom_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(11), top_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(11), back_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(11), front_side));

  assert!( rmesh3x4x5.is_boundary_side(FENum(12), left_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), right_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(12), bottom_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), top_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), back_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), front_side));

  assert!( rmesh3x4x5.is_boundary_side(FENum(4*12), left_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(4*12), right_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(4*12), bottom_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(4*12), top_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(4*12), back_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(4*12), front_side));

  assert!(!rmesh3x4x5.is_boundary_side(FENum(5*12-1), left_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(5*12-1), right_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(5*12-1), bottom_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(5*12-1), top_side));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(5*12-1), back_side));
  assert!( rmesh3x4x5.is_boundary_side(FENum(5*12-1), front_side));
}

fn test_3x4x5x6_boundary_side_determ() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));
  let early_side = lesser_side_face_perp_to_axis(Dim(3));
  let late_side = greater_side_face_perp_to_axis(Dim(3));

  assert!( rmesh3x4x5x6.is_boundary_side(FENum(0), left_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(0), right_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(0), bottom_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(0), top_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(0), back_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(0), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(0), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(0), late_side));

  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(2), left_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(2), right_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(2), bottom_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(2), top_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(2), back_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(2), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(2), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(2), late_side));

  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3), left_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3), right_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3), bottom_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3), top_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3), back_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3), late_side));

  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(11), left_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(11), right_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(11), bottom_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(11), top_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(11), back_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(11), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(11), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(11), late_side));

  assert!( rmesh3x4x5x6.is_boundary_side(FENum(12), left_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(12), right_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(12), bottom_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(12), top_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(12), back_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(12), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(12), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(12), late_side));

  assert!( rmesh3x4x5x6.is_boundary_side(FENum(4*12), left_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(4*12), right_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(4*12), bottom_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(4*12), top_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(4*12), back_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(4*12), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(4*12), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(4*12), late_side));

  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), left_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), right_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), bottom_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), top_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), back_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), front_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12-1), late_side));

  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12), left_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12), right_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12), bottom_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12), top_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(5*12), back_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12), front_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12), early_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(5*12), late_side));
  
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), left_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), right_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), bottom_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), top_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), back_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), front_side));
  assert!(!rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), early_side));
  assert!( rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6-1), late_side));
}

#[test]
#[should_fail]
fn test_3x4_bad_is_boundary_side_fenum() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  rmesh3x4.is_boundary_side(FENum(12), top_side);
}

#[test]
#[should_fail]
fn test_3x4x5_bad_is_boundary_side_fenum() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let front_side = greater_side_face_perp_to_axis(Dim(2));
  rmesh3x4x5.is_boundary_side(FENum(5*12), front_side);
}

#[test]
#[should_fail]
fn test_3x4x5x6_bad_is_boundary_side_fenum() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3.],
                                                     ~[2f64, 3., 4.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let front_side = greater_side_face_perp_to_axis(Dim(2));
  rmesh3x4x5x6.is_boundary_side(FENum(3*4*5*6), front_side);
}

// Test the non-boundary sides perpendicular to axis 0 for 2d mesh.
#[test]
fn test_3x4_nonboundary_side_coords_axis0() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));

  // The mesh for non-boundary sides perpendicular to axis 0 has dimensions 2 x 4.

  // first side in first row (axis 0)
  let sgeom_0 = rmesh3x4.nb_side_geom(NBSideNum(0));
  assert_eq!(sgeom_0.perp_axis, Dim(0));
  assert_eq!(&sgeom_0.mesh_coords, &~[MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(0));
  assert_eq!(rmesh3x4.fe_inclusions_of_nb_side(NBSideNum(0)),
             NBSideInclusions { nb_side_num: NBSideNum(0),
                                fe1: FENum(0), sideface_in_fe1: right_side,
                                fe2: FENum(1), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(0), right_side), NBSideNum(0));
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(1), left_side),  NBSideNum(0));

  // second side in first row
  let sgeom_1 = rmesh3x4.nb_side_geom(NBSideNum(1));
  assert_eq!(sgeom_1.perp_axis, Dim(0));
  assert_eq!(&sgeom_1.mesh_coords, &~[MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(0)], Dim(0)), NBSideNum(1));
  assert_eq!(rmesh3x4.fe_inclusions_of_nb_side(NBSideNum(1)),
             NBSideInclusions { nb_side_num: NBSideNum(1),
                                fe1: FENum(1), sideface_in_fe1: right_side,
                                fe2: FENum(2), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(1), right_side), NBSideNum(1));
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(2), left_side),  NBSideNum(1));

  // first side in second row
  let sgeom_2 = rmesh3x4.nb_side_geom(NBSideNum(2));
  assert_eq!(sgeom_2.perp_axis, Dim(0));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1)], Dim(0)), NBSideNum(2));
  assert_eq!(rmesh3x4.fe_inclusions_of_nb_side(NBSideNum(2)),
             NBSideInclusions { nb_side_num: NBSideNum(2),
                                fe1: FENum(3), sideface_in_fe1: right_side,
                                fe2: FENum(4), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(3), right_side), NBSideNum(2));
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(4), left_side),  NBSideNum(2));

  // last side
  let sgeom_7 = rmesh3x4.nb_side_geom(NBSideNum(7));
  assert_eq!(sgeom_7.perp_axis, Dim(0));
  assert_eq!(&sgeom_7.mesh_coords, &~[MeshCoord(1), MeshCoord(3)]);
  assert_eq!(rmesh3x4.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3)], Dim(0)), NBSideNum(7));
  assert_eq!(rmesh3x4.fe_inclusions_of_nb_side(NBSideNum(7)),
             NBSideInclusions { nb_side_num: NBSideNum(7),
                                fe1: FENum(10), sideface_in_fe1: right_side,
                                fe2: FENum(11), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(10), right_side), NBSideNum(7));
  assert_eq!(rmesh3x4.nb_side_num_for_fe_side(FENum(11), left_side),  NBSideNum(7));
}


// Test the non-boundary sides perpendicular to axis 0 for 3d mesh.
#[test]
fn test_3x4x5_nonboundary_side_coords_axis0() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));

  // The mesh for non-boundary sides perpendicular to axis 0 has dimensions 2 x 4 x 5.

  // first side in first row (axis 0)
  let sgeom_0 = rmesh3x4x5.nb_side_geom(NBSideNum(0));
  assert_eq!(sgeom_0.perp_axis, Dim(0));
  assert_eq!(&sgeom_0.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(0));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(0)),
             NBSideInclusions { nb_side_num: NBSideNum(0),
                                fe1: FENum(0), sideface_in_fe1: right_side,
                                fe2: FENum(1), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(0), right_side), NBSideNum(0));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(1), left_side),  NBSideNum(0));

  // second side in first row
  let sgeom_1 = rmesh3x4x5.nb_side_geom(NBSideNum(1));
  assert_eq!(sgeom_1.perp_axis, Dim(0));
  assert_eq!(&sgeom_1.mesh_coords, &~[MeshCoord(1), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(1));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(1)),
             NBSideInclusions { nb_side_num: NBSideNum(1),
                                fe1: FENum(1), sideface_in_fe1: right_side,
                                fe2: FENum(2), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(1), right_side), NBSideNum(1));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(2), left_side),  NBSideNum(1));

  // first side in second row
  let sgeom_2 = rmesh3x4x5.nb_side_geom(NBSideNum(2));
  assert_eq!(sgeom_2.perp_axis, Dim(0));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(0)), NBSideNum(2));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(2)),
             NBSideInclusions { nb_side_num: NBSideNum(2),
                                fe1: FENum(3), sideface_in_fe1: right_side,
                                fe2: FENum(4), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3), right_side), NBSideNum(2));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(4), left_side),  NBSideNum(2));

  // last side in first 2-stack
  let sgeom_7 = rmesh3x4x5.nb_side_geom(NBSideNum(7));
  assert_eq!(sgeom_7.perp_axis, Dim(0));
  assert_eq!(&sgeom_7.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(0)], Dim(0)), NBSideNum(7));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(7)),
             NBSideInclusions { nb_side_num: NBSideNum(7),
                                fe1: FENum(10), sideface_in_fe1: right_side,
                                fe2: FENum(11), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(10), right_side), NBSideNum(7));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(11), left_side),  NBSideNum(7));

  // first side of second 2-stack
  let sgeom_8 = rmesh3x4x5.nb_side_geom(NBSideNum(8));
  assert_eq!(sgeom_8.perp_axis, Dim(0));
  assert_eq!(&sgeom_8.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(0)), NBSideNum(8));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(8)),
             NBSideInclusions { nb_side_num: NBSideNum(8),
                                fe1: FENum(12), sideface_in_fe1: right_side,
                                fe2: FENum(13), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), right_side), NBSideNum(8));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(13), left_side),  NBSideNum(8));

  // first side of last 2-stack
  let sgeom_32 = rmesh3x4x5.nb_side_geom(NBSideNum(32));
  assert_eq!(sgeom_32.perp_axis, Dim(0));
  assert_eq!(&sgeom_32.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4)], Dim(0)), NBSideNum(32));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(32)),
             NBSideInclusions { nb_side_num: NBSideNum(32),
                                fe1: FENum(48), sideface_in_fe1: right_side,
                                fe2: FENum(49), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), right_side), NBSideNum(32));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(49), left_side),  NBSideNum(32));

  // last side of last 2-stack
  let sgeom_39 = rmesh3x4x5.nb_side_geom(NBSideNum(39));
  assert_eq!(sgeom_39.perp_axis, Dim(0));
  assert_eq!(&sgeom_39.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(4)], Dim(0)), NBSideNum(39));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(39)),
             NBSideInclusions { nb_side_num: NBSideNum(39),
                                fe1: FENum(58), sideface_in_fe1: right_side,
                                fe2: FENum(59), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(58), right_side), NBSideNum(39));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(59), left_side),  NBSideNum(39));
}

// Test the non-boundary sides perpendicular to axis 0 for 4d mesh.
#[test]
fn test_3x4x5x6_nonboundary_side_coords_axis0() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));

  // The mesh for non-boundary sides perpendicular to axis 0 has dimensions 2 x 4 x 5 x 6.

  // first side in first row (axis 0)
  let sgeom_0 = rmesh3x4x5x6.nb_side_geom(NBSideNum(0));
  assert_eq!(sgeom_0.perp_axis, Dim(0));
  assert_eq!(&sgeom_0.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(0));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(0)),
             NBSideInclusions { nb_side_num: NBSideNum(0),
                                fe1: FENum(0), sideface_in_fe1: right_side,
                                fe2: FENum(1), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(0), right_side), NBSideNum(0));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(1), left_side),  NBSideNum(0));

  // second side in first row
  let sgeom_1 = rmesh3x4x5x6.nb_side_geom(NBSideNum(1));
  assert_eq!(sgeom_1.perp_axis, Dim(0));
  assert_eq!(&sgeom_1.mesh_coords, &~[MeshCoord(1), MeshCoord(0), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(0), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(1));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(1)),
             NBSideInclusions { nb_side_num: NBSideNum(1),
                                fe1: FENum(1), sideface_in_fe1: right_side,
                                fe2: FENum(2), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(1), right_side), NBSideNum(1));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(2), left_side),  NBSideNum(1));

  // first side in second row
  let sgeom_2 = rmesh3x4x5x6.nb_side_geom(NBSideNum(2));
  assert_eq!(sgeom_2.perp_axis, Dim(0));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(2));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(2)),
             NBSideInclusions { nb_side_num: NBSideNum(2),
                                fe1: FENum(3), sideface_in_fe1: right_side,
                                fe2: FENum(4), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(3), right_side), NBSideNum(2));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(4), left_side),  NBSideNum(2));

  // last side in first 2-stack
  let sgeom_7 = rmesh3x4x5x6.nb_side_geom(NBSideNum(7));
  assert_eq!(sgeom_7.perp_axis, Dim(0));
  assert_eq!(&sgeom_7.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(7));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(7)),
             NBSideInclusions { nb_side_num: NBSideNum(7),
                                fe1: FENum(10), sideface_in_fe1: right_side,
                                fe2: FENum(11), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(10), right_side), NBSideNum(7));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(11), left_side),  NBSideNum(7));

  // first side of second 2-stack
  let sgeom_8 = rmesh3x4x5x6.nb_side_geom(NBSideNum(8));
  assert_eq!(sgeom_8.perp_axis, Dim(0));
  assert_eq!(&sgeom_8.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(0)), NBSideNum(8));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(8)),
             NBSideInclusions { nb_side_num: NBSideNum(8),
                                fe1: FENum(12), sideface_in_fe1: right_side,
                                fe2: FENum(13), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(12), right_side), NBSideNum(8));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(13), left_side),  NBSideNum(8));

  // first side of last 2-stack
  let sgeom_32 = rmesh3x4x5x6.nb_side_geom(NBSideNum(32));
  assert_eq!(sgeom_32.perp_axis, Dim(0));
  assert_eq!(&sgeom_32.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4), MeshCoord(0)], Dim(0)), NBSideNum(32));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(32)),
             NBSideInclusions { nb_side_num: NBSideNum(32),
                                fe1: FENum(48), sideface_in_fe1: right_side,
                                fe2: FENum(49), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(48), right_side), NBSideNum(32));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(49), left_side),  NBSideNum(32));

  // last side of last 2-stack
  let sgeom_39 = rmesh3x4x5x6.nb_side_geom(NBSideNum(39));
  assert_eq!(sgeom_39.perp_axis, Dim(0));
  assert_eq!(&sgeom_39.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(4), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(4), MeshCoord(0)], Dim(0)), NBSideNum(39));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(39)),
             NBSideInclusions { nb_side_num: NBSideNum(39),
                                fe1: FENum(58), sideface_in_fe1: right_side,
                                fe2: FENum(59), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(58), right_side), NBSideNum(39));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(59), left_side),  NBSideNum(39));

  // first side of second 3-stack
  let sgeom_40 = rmesh3x4x5x6.nb_side_geom(NBSideNum(40));
  assert_eq!(sgeom_40.perp_axis, Dim(0));
  assert_eq!(&sgeom_40.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(0)), NBSideNum(40));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(40)),
             NBSideInclusions { nb_side_num: NBSideNum(40),
                                fe1: FENum(60), sideface_in_fe1: right_side,
                                fe2: FENum(61), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(60), right_side), NBSideNum(40));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(61), left_side),  NBSideNum(40));
  
  // second side of second 3-stack
  let sgeom_41 = rmesh3x4x5x6.nb_side_geom(NBSideNum(41));
  assert_eq!(sgeom_41.perp_axis, Dim(0));
  assert_eq!(&sgeom_41.mesh_coords, &~[MeshCoord(1), MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(0)), NBSideNum(41));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(41)),
             NBSideInclusions { nb_side_num: NBSideNum(41),
                                fe1: FENum(61), sideface_in_fe1: right_side,
                                fe2: FENum(62), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(61), right_side), NBSideNum(41));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(62), left_side),  NBSideNum(41));
  
  // last side
  let sgeom_last = rmesh3x4x5x6.nb_side_geom(NBSideNum(2*4*5*6-1));
  assert_eq!(sgeom_last.perp_axis, Dim(0));
  assert_eq!(&sgeom_last.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  assert_eq!(rmesh3x4x5x6.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(4), MeshCoord(5)], Dim(0)), NBSideNum(2*4*5*6-1));
  assert_eq!(rmesh3x4x5x6.fe_inclusions_of_nb_side(NBSideNum(2*4*5*6-1)),
             NBSideInclusions { nb_side_num: NBSideNum(2*4*5*6-1),
                                fe1: FENum(3*4*5*6-2), sideface_in_fe1: right_side,
                                fe2: FENum(3*4*5*6-1), sideface_in_fe2: left_side });
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(3*4*5*6-2), right_side), NBSideNum(2*4*5*6-1));
  assert_eq!(rmesh3x4x5x6.nb_side_num_for_fe_side(FENum(3*4*5*6-1), left_side),  NBSideNum(2*4*5*6-1));
}

// Test the non-boundary sides perpendicular to axis 1 (3d mesh).
#[test]
fn test_3x4x5_nonboundary_side_coords_axis1() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));

  let first_axis1 = 2*4*5;
  assert_eq!(first_axis1, *rmesh3x4x5.first_nb_side_nums_by_perp_axis[1]);

  // The mesh for non-boundary sides perpendicular to axis 1 has dimensions 3 x 3 x 5.

  // first side in first row (axis 1)
  let sgeom_0 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+0));
  assert_eq!(sgeom_0.perp_axis, Dim(1));
  assert_eq!(&sgeom_0.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+0));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+0)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+0),
                                fe1: FENum(0), sideface_in_fe1: top_side,
                                fe2: FENum(3), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(0), top_side), NBSideNum(first_axis1+0));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3), bottom_side),  NBSideNum(first_axis1+0));

  // last side in first row
  let sgeom_2 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+2));
  assert_eq!(sgeom_2.perp_axis, Dim(1));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(2), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(0), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+2));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+2)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+2),
                                fe1: FENum(2), sideface_in_fe1: top_side,
                                fe2: FENum(5), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(2), top_side), NBSideNum(first_axis1+2));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(5), bottom_side),  NBSideNum(first_axis1+2));

  // first side in second row
  let sgeom_3 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+3));
  assert_eq!(sgeom_3.perp_axis, Dim(1));
  assert_eq!(&sgeom_3.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+3));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+3)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+3),
                                fe1: FENum(3), sideface_in_fe1: top_side,
                                fe2: FENum(6), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3), top_side), NBSideNum(first_axis1+3));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(6), bottom_side),  NBSideNum(first_axis1+3));

  // last side in first stack
  let sgeom_8 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+8));
  assert_eq!(sgeom_8.perp_axis, Dim(1));
  assert_eq!(&sgeom_8.mesh_coords, &~[MeshCoord(2), MeshCoord(2), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(2), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+8));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+8)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+8),
                                fe1: FENum(8), sideface_in_fe1: top_side,
                                fe2: FENum(11), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(8), top_side), NBSideNum(first_axis1+8));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(11), bottom_side),  NBSideNum(first_axis1+8));

  // first side of second stack
  let sgeom_9 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+9));
  assert_eq!(sgeom_9.perp_axis, Dim(1));
  assert_eq!(&sgeom_9.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(1)), NBSideNum(first_axis1+9));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+9)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+9),
                                fe1: FENum(12), sideface_in_fe1: top_side,
                                fe2: FENum(15), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), top_side), NBSideNum(first_axis1+9));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(15), bottom_side),  NBSideNum(first_axis1+9));

  // first side of last stack
  let sgeom_36 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+36));
  assert_eq!(sgeom_36.perp_axis, Dim(1));
  assert_eq!(&sgeom_36.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4)], Dim(1)), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+36)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+36),
                                fe1: FENum(48), sideface_in_fe1: top_side,
                                fe2: FENum(51), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), top_side), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(51), bottom_side),  NBSideNum(first_axis1+36));

  // last side of last stack
  let sgeom_36 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+36));
  assert_eq!(sgeom_36.perp_axis, Dim(1));
  assert_eq!(&sgeom_36.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4)], Dim(1)), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+36)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+36),
                                fe1: FENum(48), sideface_in_fe1: top_side,
                                fe2: FENum(51), sideface_in_fe2: bottom_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), top_side), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(51), bottom_side),  NBSideNum(first_axis1+36));

}

// Test the non-boundary sides perpendicular to axis 2 (3d mesh).
#[test]
fn test_3x4x5_nonboundary_side_coords_axis2() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let front_side = lesser_side_face_perp_to_axis(Dim(2));
  let back_side = greater_side_face_perp_to_axis(Dim(2));

  let first_axis2 = 2*4*5 + 3*3*5;
  assert_eq!(first_axis2, *rmesh3x4x5.first_nb_side_nums_by_perp_axis[2]);

  // The mesh for non-boundary sides perpendicular to axis 2 has dimensions 3 x 4 x 4.

  // first side in first row (axis 2)
  let sgeom_0 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+0));
  assert_eq!(sgeom_0.perp_axis, Dim(2));
  assert_eq!(&sgeom_0.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(0)], Dim(2)), NBSideNum(first_axis2+0));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+0)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+0),
                                fe1: FENum(0),  sideface_in_fe1: back_side,
                                fe2: FENum(12), sideface_in_fe2: front_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(0),  back_side), NBSideNum(first_axis2+0));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), front_side),  NBSideNum(first_axis2+0));

  // third side in first row
  let sgeom_2 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+2));
  assert_eq!(sgeom_2.perp_axis, Dim(2));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(2), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(0), MeshCoord(0)], Dim(2)), NBSideNum(first_axis2+2));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+2)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+2),
                                fe1: FENum(2),  sideface_in_fe1: back_side,
                                fe2: FENum(14), sideface_in_fe2: front_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(2),  back_side), NBSideNum(first_axis2+2));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(14), front_side),  NBSideNum(first_axis2+2));


  // second row
  let sgeom_3 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+3));
  assert_eq!(sgeom_3.perp_axis, Dim(2));
  assert_eq!(&sgeom_3.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(2)), NBSideNum(first_axis2+3));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+3)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+3),
                                fe1: FENum(3),  sideface_in_fe1: back_side,
                                fe2: FENum(15), sideface_in_fe2: front_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3),  back_side), NBSideNum(first_axis2+3));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(15), front_side),  NBSideNum(first_axis2+3));
  
  // first side of second stack
  let sgeom_12 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+12));
  assert_eq!(sgeom_12.perp_axis, Dim(2));
  assert_eq!(&sgeom_12.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(2)), NBSideNum(first_axis2+12));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+12)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+12),
                                fe1: FENum(12), sideface_in_fe1: back_side,
                                fe2: FENum(24), sideface_in_fe2: front_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), back_side), NBSideNum(first_axis2+12));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(24), front_side),  NBSideNum(first_axis2+12));

  // first side of last stack
  let sgeom_36 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+36));
  assert_eq!(sgeom_36.perp_axis, Dim(2));
  assert_eq!(&sgeom_36.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(3)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(3)], Dim(2)), NBSideNum(first_axis2+36));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+36)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+36),
                                fe1: FENum(36), sideface_in_fe1: back_side,
                                fe2: FENum(48), sideface_in_fe2: front_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(36), back_side), NBSideNum(first_axis2+36));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), front_side),  NBSideNum(first_axis2+36));

  // last side of last stack
  let sgeom_47 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+47));
  assert_eq!(sgeom_47.perp_axis, Dim(2));
  assert_eq!(&sgeom_47.mesh_coords, &~[MeshCoord(2), MeshCoord(3), MeshCoord(3)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(3), MeshCoord(3)], Dim(2)), NBSideNum(first_axis2+47));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+47)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+47),
                                fe1: FENum(47), sideface_in_fe1: back_side,
                                fe2: FENum(59), sideface_in_fe2: front_side });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(47), back_side), NBSideNum(first_axis2+47));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(59), front_side), NBSideNum(first_axis2+47));

  // Above was the last non-boundary side.  Check that the total number is correct.
  assert_eq!(rmesh3x4x5.num_boundary_sides(), 2*20 + 2*15 + 2*12);
}

// Attempt to access beyond the last non-boundary side should fail (3d mesh).
#[test]
#[should_fail]
fn test_3x4x5_nonboundary_bad_side_coords_axis2() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let first_axis2 = 2*4*5 + 3*3*5;
  rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+48));
}


#[test]
fn test_3x4_conv_fe_coords_to_fenum() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  assert_eq!(rmesh3x4.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0)]), FENum(0));
  assert_eq!(rmesh3x4.fe_with_mesh_coords([MeshCoord(2), MeshCoord(0)]), FENum(2));
  assert_eq!(rmesh3x4.fe_with_mesh_coords([MeshCoord(0), MeshCoord(1)]), FENum(3));
  assert_eq!(rmesh3x4.fe_with_mesh_coords([MeshCoord(1), MeshCoord(1)]), FENum(4));
}

#[test]
fn test_3x4x5_conv_fe_coords_to_fenum() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(0)]), FENum(0));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(2), MeshCoord(0), MeshCoord(0)]), FENum(2));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(0), MeshCoord(1), MeshCoord(0)]), FENum(3));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(1), MeshCoord(1), MeshCoord(0)]), FENum(4));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(1)]), FENum(12));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(2), MeshCoord(3), MeshCoord(4)]), FENum(59));
}

#[test]
fn test_3x4x5x6_conv_fe_coords_to_fenum() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                   ~[2f64, 3., 4., 5.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(0), MeshCoord(0)]), FENum(0));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(2), MeshCoord(0), MeshCoord(0), MeshCoord(0)]), FENum(2));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(0), MeshCoord(1), MeshCoord(0), MeshCoord(0)]), FENum(3));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(1), MeshCoord(1), MeshCoord(0), MeshCoord(0)]), FENum(4));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(1), MeshCoord(0)]), FENum(12));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(2), MeshCoord(3), MeshCoord(4), MeshCoord(0)]), FENum(59));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(0), MeshCoord(1)]), FENum(60));
  assert_eq!(rmesh3x4x5x6.fe_with_mesh_coords([MeshCoord(2), MeshCoord(3), MeshCoord(4), MeshCoord(5)]), FENum(3*4*5*6-1));
}

// Test integration of monomials through RectIntegrable trait.

#[test]
fn test_intg_mon1d() {
  let one = Mon1d { exps: [Deg(0)] };
  let x = Mon1d { exps: [Deg(1)] };
  let x2 = Mon1d { exps: [Deg(2)] };
  assert_eq!(one.integral_over_rect_at_origin([2.]), 2.);
  assert_eq!(x.integral_over_rect_at_origin([2.]), 2.);
  assert_eq!(x2.integral_over_rect_at_origin([2.]), 8./3.);
}

#[test]
fn test_intg_mon2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x1y2 = Mon2d { exps: [Deg(1), Deg(2)] };
  let x3y1 = Mon2d { exps: [Deg(3), Deg(1)] };
  assert_eq!(one.integral_over_rect_at_origin([2.,3.]), 6.);
  assert_eq!(x1y2.integral_over_rect_at_origin([2.,3.]), 18.);
  assert_eq!(x3y1.integral_over_rect_at_origin([2.,3.]), 4.*9./2.);
}

#[test]
fn test_intg_mon3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x1y2z3 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  let x3y1z2 = Mon3d { exps: [Deg(3), Deg(1), Deg(2)] };
  assert_eq!(one.integral_over_rect_at_origin([2.,3.,4.]), 24.);
  assert_eq!(x1y2z3.integral_over_rect_at_origin([2.,3.,4.]), 18.*64.);
  assert_eq!(x3y1z2.integral_over_rect_at_origin([2.,3.,4.]), 4.*(9./2.)*64./3.);
}

#[test]
fn test_intg_mon4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  let x1y2z3t4 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  let x3y1z2t4 = Mon4d { exps: [Deg(3), Deg(1), Deg(2), Deg(4)] };
  assert_eq!(one.integral_over_rect_at_origin([2.,3.,4.,5.]), 24.*5.);
  assert_eq!(x1y2z3t4.integral_over_rect_at_origin([2.,3.,4.,5.]), 18.*64.*625.);
  assert_eq!(x3y1z2t4.integral_over_rect_at_origin([2.,3.,4.,5.]), 4.*(9./2.)*(64./3.)*625.);
}

#[test]
fn test_intg_const_facerel_mons_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let one = Mon2d { exps: [Deg(0), Deg(0)] };

  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));

  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_int(one, OShape(0)),
                1./3. * 1./4.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(one, OShape(0), left_side),
                1./4.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(one, OShape(0), right_side),
                1./4.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(one, OShape(0), bottom_side),
                1./3.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(one, OShape(0), top_side),
                1./3.);
}

#[test]
fn test_intg_const_facerel_mons_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };

  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_int(one, OShape(0)),
                1./3. * 1./4. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(one, OShape(0), left_side),
                1./4. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(one, OShape(0), right_side),
                1./4. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(one, OShape(0), bottom_side),
                1./3. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(one, OShape(0), top_side),
                1./3. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(one, OShape(0), back_side),
                1./3. * 1./4.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(one, OShape(0), front_side),
                1./3. * 1./4.);
}

#[test]
fn test_intg_const_facerel_mons_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };

  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));
  let early_side = lesser_side_face_perp_to_axis(Dim(3));
  let late_side = greater_side_face_perp_to_axis(Dim(3));
  
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_int(one, OShape(0)),
                1./3. * 1./4. * 1./5. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), left_side),
                1./4. * 1./5.* 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), right_side),
                1./4. * 1./5. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), bottom_side),
                1./3. * 1./5. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), top_side),
                1./3. * 1./5. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), back_side),
                1./3. * 1./4. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), front_side),
                1./3. * 1./4. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), early_side),
                1./3. * 1./4. * 1./5.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(one, OShape(0), late_side),
                1./3. * 1./4. * 1./5.);
}

#[test]
fn test_intg_global_fn_on_fe_int_2d() -> () {
  let mins = ~[1f64, 2.];
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(mins.clone(),
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let x3_y4 = |x:&[R]| -> R {
    pow(x[0]-mins[0], 3) * pow(x[1]-mins[1], 4)
  };
  
  assert_approx(rmesh3x4.intg_global_fn_on_fe_int(x3_y4, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
}

#[test]
fn test_intg_global_fn_on_fe_int_3d() -> () {
  let mins = ~[1f64, 2., 3.];
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(mins.clone(),
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x3_y4_z = |x:&[R]| -> R {
    pow(x[0]-mins[0], 3) * pow(x[1]-mins[1], 4) * (x[2]-mins[2])
  };
  
  assert_approx(rmesh3x4x5.intg_global_fn_on_fe_int(x3_y4_z, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,2)/2.);
}

#[test]
fn test_intg_global_fn_on_fe_int_4d() -> () {
  let mins = ~[1f64, 2., 3., 4.];
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(mins.clone(),
                                                   ~[2f64, 3., 4., 5.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let x3_y4_z2_t = |x:&[R]| -> R {
    pow(x[0]-mins[0], 3) * pow(x[1]-mins[1], 4) * pow(x[2]-mins[2], 2) * (x[3]-mins[3])
  };
  
  assert_approx(rmesh3x4x5x6.intg_global_fn_on_fe_int(x3_y4_z2_t, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,3)/3. * 1./6.);
}

#[test]
fn test_intg_global_fn_x_facerel_mon_on_fe0_int_2d() -> () {
  let mins = ~[1f64, 2.];
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(mins.clone(),
                                                   ~[2f64, 3.],
                                                   ~[MeshCoord(3), MeshCoord(4)]);
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  let x2_y3 = |x: &[R]| -> R {
    pow(x[0]-mins[0], 2) * pow(x[1]-mins[1], 3)
  };

  assert_approx(rmesh3x4.intg_global_fn_x_facerel_mon_on_fe_int(x2_y3, x*y, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
}

#[test]
fn test_intg_global_fn_x_facerel_mon_on_fe0_int_3d() -> () {
  let mins = ~[1f64, 2., 3.];
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(mins.clone(),
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  let x2_y3 = |x: &[R]| -> R {
    pow(x[0]-mins[0], 2) * pow(x[1]-mins[1], 3)
  };

  assert_approx(rmesh3x4x5.intg_global_fn_x_facerel_mon_on_fe_int(x2_y3, x*y*z, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,2)/2.);
}

#[test]
fn test_intg_global_fn_x_facerel_mon_on_fe0_int_4d() -> () {
  let mins = ~[1f64, 2., 3., 4.];
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(mins.clone(),
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  let x2_y3 = |x: &[R]| -> R {
    pow(x[0]-mins[0], 2) * pow(x[1]-mins[1], 3)
  };

  assert_approx(rmesh3x4x5x6.intg_global_fn_x_facerel_mon_on_fe_int(x2_y3, x*y*z*z*t, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,3)/3. * pow(1./6.,2)/2.);
}

#[test]
fn test_intg_global_fn_x_facerel_mon_on_fe4_int_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  let fe4 = rmesh3x4.fe_interior_origin(FENum(4));
  let x2_y3 = |x: &[R]| -> R {
    pow(x[0]-fe4[0], 2) * pow(x[1]-fe4[1], 3)
  };

  assert_approx(rmesh3x4.intg_global_fn_x_facerel_mon_on_fe_int(x2_y3, x*y, FENum(4)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
}

#[test]
fn test_intg_global_fn_x_facerel_mon_on_fe16_int_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  let fe16 = rmesh3x4x5.fe_interior_origin(FENum(16));
  let x2_y3 = |x: &[R]| -> R {
    pow(x[0]-fe16[0], 2) * pow(x[1]-fe16[1], 3)
  };

  assert_approx(rmesh3x4x5.intg_global_fn_x_facerel_mon_on_fe_int(x2_y3, x*y*z, FENum(16)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,2)/2.);
}

#[test]
fn test_intg_global_fn_x_facerel_mon_on_fe16_int_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  let fe16 = rmesh3x4x5x6.fe_interior_origin(FENum(16));
  let x2_y3 = |x: &[R]| -> R {
    pow(x[0]-fe16[0], 2) * pow(x[1]-fe16[1], 3)
  };

  assert_approx(rmesh3x4x5x6.intg_global_fn_x_facerel_mon_on_fe_int(x2_y3, y*z*t, FENum(16)),
                pow(1./3.,3)/3. * pow(1./4.,5)/5. * pow(1./5.,2)/2. * pow(1./6.,2)/2.);
}

#[test]
fn test_intg_facerel_poly_on_oshape_int_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_approx(rmesh3x4.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y), (3.,x*x*x*y*y*y*y)]), OShape(0)),
                2.*(1./3. * pow(1./4.,2)/2.) + 3.*(pow(1./3.,4)/4. * pow(1./4.,5)/5.));

  assert_approx(rmesh3x4.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y)]), OShape(0)),
                1./3. * pow(1./4.,2));

  assert_approx(rmesh3x4.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y), (3.,x*x)]), OShape(0)),
                2. * 1./3. * pow(1./4.,2)/2.  + 3. * pow(1./3.,3)/3. * 1./4.);
}

#[test]
fn test_intg_facerel_poly_on_oshape_int_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y*z*z), (3.,x*x*x*y*y*y*y*z)]), OShape(0)),
                2.*(1./3. * pow(1./4.,2)/2. * pow(1./5.,3)/3.) + 3.*(pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,2)/2.));

  assert_approx(rmesh3x4x5.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y)]), OShape(0)),
                1./3. * pow(1./4.,2) * 1./5.);

  assert_approx(rmesh3x4x5.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y), (3.,x*x)]), OShape(0)),
                2. * 1./3. * pow(1./4.,2)/2. * 1./5. + 3. * pow(1./3.,3)/3. * 1./4. * 1./5.);
}

#[test]
fn test_intg_facerel_poly_on_oshape_int_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5x6.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y*z*z*t), (3.,x*x*x*y*y*y*y*z)]), OShape(0)),
                2.*(1./3. * pow(1./4.,2)/2. * pow(1./5.,3)/3. * pow(1./6.,2)/2.) + 3.*(pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,2)/2. * 1./6.));

  assert_approx(rmesh3x4x5x6.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y)]), OShape(0)),
                1./3. * pow(1./4.,2) * 1./5. * 1./6.);

  assert_approx(rmesh3x4x5x6.intg_facerel_poly_on_oshape_int(&poly(~[(2.,y), (3.,x*x*t)]), OShape(0)),
                2. * 1./3. * pow(1./4.,2)/2. * 1./5. * 1./6. + 3. * pow(1./3.,3)/3. * 1./4. * 1./5. * pow(1./6.,2)/2.);
}

#[test]
fn test_intg_facerel_poly_x_facerel_poly_on_oshape_int_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                   ~[2f64, 3.],
                                                   ~[MeshCoord(3), MeshCoord(4)]);
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  // (2y + 3xy)(x - 2y) = 2xy - 4y^2 + 3x^2y - 6xy^2
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_int(&poly(~[(2.,y), (3.,x*y)]), &poly(~[(1.,x), (-2.,y)]), OShape(0)),
                2. * pow(1./3.,2)/2. * pow(1./4.,2)/2.
              - 4. * 1./3.           * pow(1./4.,3)/3.
              + 3. * pow(1./3.,3)/3. * pow(1./4.,2)/2.
              - 6. * pow(1./3.,2)/2. * pow(1./4.,3)/3.);
}

#[test]
fn test_intg_facerel_poly_x_facerel_poly_on_oshape_int_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  // (2yz^2 + 3xyz)(x - 2y) = 2xyz^2 - 4y^2z^2 + 3x^2yz - 6xy^2z
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_int(&poly(~[(2.,y*z*z), (3.,x*y*z)]), &poly(~[(1.,x), (-2.,y)]), OShape(0)),
                2. * pow(1./3.,2)/2. * pow(1./4.,2)/2. * pow(1./5.,3)/3. 
              - 4. * 1./3.           * pow(1./4.,3)/3. * pow(1./5.,3)/3.
              + 3. * pow(1./3.,3)/3. * pow(1./4.,2)/2. * pow(1./5.,2)/2. 
              - 6. * pow(1./3.,2)/2. * pow(1./4.,3)/3. * pow(1./5.,2)/2.);
}

#[test]
fn test_intg_facerel_poly_x_facerel_poly_on_oshape_int_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };

  // (2yz^2 + 3xyz)(x - 2y) = 2xyz^2 - 4y^2z^2 + 3x^2yz - 6xy^2z
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_int(&poly(~[(2.,y*z*z), (3.,x*y*z)]), &poly(~[(1.,x), (-2.,y)]), OShape(0)),
                2. * pow(1./3.,2)/2. * pow(1./4.,2)/2. * pow(1./5.,3)/3. * 1./6.
              - 4. * 1./3.           * pow(1./4.,3)/3. * pow(1./5.,3)/3. * 1./6.
              + 3. * pow(1./3.,3)/3. * pow(1./4.,2)/2. * pow(1./5.,2)/2. * 1./6.
              - 6. * pow(1./3.,2)/2. * pow(1./4.,3)/3. * pow(1./5.,2)/2. * 1./6.);
}

fn test_intg_facerel_poly_x_facerel_poly_on_oshape_side_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2.],
                                                   ~[2f64, 3.],
                                                   ~[MeshCoord(3), MeshCoord(4)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let one_poly = ~poly(~[(1.,one)]);

  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.2,y), (3.4,x*x*x*y*y*y*y)]), one_poly, OShape(0), left_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.2,y), (3.4,x*x*x*y*y*y*y)]), one_poly, OShape(0), right_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
 
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.5,y), (0.123,y*y*y*y)]), one_poly, OShape(0), right_side),
                1.5 * pow(1./4.,2)/2. * pow(1./5.,3)/3. + 0.123 * pow(1./4.,5)/5. * pow(1./5.,2)/2.);

  // (2y - 3)(2 - y^2) = 4y - 2y^3 - 6 + 3y^2
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,y), (-3.,one)]), &poly(~[(2.,one), (-1.,y*y)]), OShape(0), right_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3.
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1.
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3.);
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,y), (-3.,one)]), &poly(~[(2.,one), (-1.,y*y)]), OShape(0), left_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3.
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1.
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3.);
  
  // (2x - 3)(2 - x^2) = 4x - 2x^3 - 6 + 3x^2
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,one), (-1.,x*x)]), OShape(0), top_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3.
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1.
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3.);
  assert_approx(rmesh3x4.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,one), (-1.,x*x)]), OShape(0), bottom_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3.
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1.
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3.);
}

fn test_intg_facerel_poly_x_facerel_poly_on_oshape_side_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_poly = ~poly(~[(1.,one)]);

  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.2,y*z*z), (3.4,x*x*x*y*y*y*y)]), one_poly, OShape(0), left_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.2,y*z*z), (3.4,x*x*x*y*y*y*y)]), one_poly, OShape(0), right_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
 
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.5,y*z*z), (0.123,y*y*y*y*z)]), one_poly, OShape(0), right_side),
                1.5 * pow(1./4.,2)/2. * pow(1./5.,3)/3. + 0.123 * pow(1./4.,5)/5. * pow(1./5.,2)/2.);

  // (2y - 3)(2z - y^2) = 4yz - 2y^3 - 6z + 3y^2
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,y), (-3.,one)]), &poly(~[(2.,z), (-1.,y*y)]), OShape(0), right_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1.);
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,y), (-3.,one)]), &poly(~[(2.,z), (-1.,y*y)]), OShape(0), left_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1.);
  
  // (2x - 3)(2z - x^2) = 4xz - 2x^3 - 6z + 3x^2
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,z), (-1.,x*x)]), OShape(0), top_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1.);
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,z), (-1.,x*x)]), OShape(0), bottom_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1.);

  // (2x - 3)(2y - x^2) = 4xy - 2x^3 - 6y + 3x^2
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,y), (-1.,x*x)]), OShape(0), front_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,1)/1. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1.);
  assert_approx(rmesh3x4x5.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,y), (-1.,x*x)]), OShape(0), back_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1.
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,1)/1. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1.);
}

fn test_intg_facerel_poly_x_facerel_poly_on_oshape_side_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));
  let early_side = lesser_side_face_perp_to_axis(Dim(3));
  let late_side = greater_side_face_perp_to_axis(Dim(3));

  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one_poly = ~poly(~[(1.,one)]);

  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.2,y*z*z*t), (3.4,x*x*x*y*y*y*y)]), one_poly, OShape(0), left_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3. * pow(1./6.,2)/2.);
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.2,y*z*z*t), (3.4,x*x*x*y*y*y*y)]), one_poly, OShape(0), right_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3. * pow(1./6.,2)/2.);
 
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(1.5,y*z*z*t), (0.123,y*y*y*y*z)]), one_poly, OShape(0), right_side),
                1.5 * pow(1./4.,2)/2. * pow(1./5.,3)/3. * pow(1./6.,2)/2. + 0.123 * pow(1./4.,5)/5. * pow(1./5.,2)/2. * 1./6.);

  // (2y - 3)(2z - y^2) = 4yz - 2y^3 - 6z + 3y^2
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,y), (-3.,one)]), &poly(~[(2.,z), (-1.,y*y)]), OShape(0), right_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1. 
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1. 
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. * pow(1./6.,1)/1. 
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,y), (-3.,one)]), &poly(~[(2.,z), (-1.,y*y)]), OShape(0), left_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
  
  // (2x - 3)(2z - x^2) = 4xz - 2x^3 - 6z + 3x^2
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,z), (-1.,x*x)]), OShape(0), top_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.   
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,z), (-1.,x*x)]), OShape(0), bottom_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,2)/2. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);

  // (2x - 3)(2y - x^2) = 4xy - 2x^3 - 6y + 3x^2
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,y), (-1.,x*x)]), OShape(0), front_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.   
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,y), (-1.,x*x)]), OShape(0), back_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.   
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
  
  // (2x - 3)(2y - x^2) = 4xy - 2x^3 - 6y + 3x^2
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,y), (-1.,x*x)]), OShape(0), early_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.   
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
  assert_approx(rmesh3x4x5x6.intg_facerel_poly_x_facerel_poly_on_oshape_side(&poly(~[(2.,x), (-3.,one)]), &poly(~[(2.,y), (-1.,x*x)]), OShape(0), late_side),
                4. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.   
              - 2. * pow(1./3.,1)/1. * pow(1./4.,4)/4. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              - 6. * pow(1./3.,1)/1. * pow(1./4.,1)/1. * pow(1./5.,1)/1. * pow(1./6.,1)/1.  
              + 3. * pow(1./3.,1)/1. * pow(1./4.,3)/3. * pow(1./5.,1)/1. * pow(1./6.,1)/1.);
}

#[test]
fn test_intg_facerel_mon_on_oshape_int_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_int(x*y, OShape(0)),
                pow(1./3.,2)/2. * pow(1./4.,2)/2.);
}


#[test]
fn test_intg_facerel_mon_on_oshape_int_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_int(x*y*z*z, OShape(0)),
                pow(1./3.,2)/2. * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
}

#[test]
fn test_intg_facerel_mon_on_oshape_int_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_int(x*y*z*z*t, OShape(0)),
                pow(1./3.,2)/2. * pow(1./4.,2)/2. * pow(1./5.,3)/3. * pow(1./6.,2)/2.);
}

#[test]
fn test_intg_facerel_mon_on_oshape_side_2d() -> () {
  let rmesh3x4: ~RectMesh<Mon2d> = RectMesh::new(~[1f64, 2.],
                                                 ~[2f64, 3.],
                                                 ~[MeshCoord(3), MeshCoord(4)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));

  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(y, OShape(0), left_side),
                pow(1./4.,2)/2.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(y, OShape(0), right_side),
                pow(1./4.,2)/2.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(y, OShape(0), bottom_side),
                0.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(y, OShape(0), top_side),
                0.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(x, OShape(0), bottom_side),
                pow(1./3.,2)/2.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(y*y*y*y, OShape(0), left_side),
                pow(1./4.,5)/5.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(y*y*y*y, OShape(0), right_side),
                pow(1./4.,5)/5.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(x*x*x, OShape(0), bottom_side),
                pow(1./3.,4)/4.);
  assert_approx(rmesh3x4.intg_facerel_mon_on_oshape_side(x*x*x, OShape(0), top_side),
                pow(1./3.,4)/4.);
}


#[test]
fn test_intg_facerel_mon_on_oshape_side_3d() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), left_side),
                pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), right_side),
                pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), bottom_side),
                0.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), top_side),
                0.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*z*z, OShape(0), bottom_side),
                pow(1./3.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*x*y, OShape(0), back_side),
                pow(1./3.,3)/3. * pow(1./4.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*x*y, OShape(0), front_side),
                pow(1./3.,3)/3. * pow(1./4.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(y*y*y*y*z, OShape(0), left_side),
                pow(1./4.,5)/5. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(y*y*y*y*z, OShape(0), right_side),
                pow(1./4.,5)/5. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*x*x*z, OShape(0), bottom_side),
                pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*x*x*z, OShape(0), top_side),
                pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*x*x*y*y*y*y, OShape(0), back_side),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_on_oshape_side(x*x*x*y*y*y*y, OShape(0), front_side),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
}

#[test]
fn test_intg_facerel_mon_on_oshape_side_4d() -> () {
  let rmesh3x4x5x6: ~RectMesh<Mon4d> = RectMesh::new(~[1f64, 2., 3., 4.],
                                                     ~[2f64, 3., 4., 5.],
                                                     ~[MeshCoord(3), MeshCoord(4), MeshCoord(5), MeshCoord(6)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));
  let early_side = lesser_side_face_perp_to_axis(Dim(3));
  let late_side = greater_side_face_perp_to_axis(Dim(3));
  
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(y*z*z*t, OShape(0), left_side),
                pow(1./4.,2)/2. * pow(1./5.,3)/3. * pow(1./6.,2)/2.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), right_side),
                pow(1./4.,2)/2. * pow(1./5.,3)/3. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), bottom_side),
                0.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(y*z*z, OShape(0), top_side),
                0.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*z*z, OShape(0), bottom_side),
                pow(1./3.,2)/2. * pow(1./5.,3)/3. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*y, OShape(0), back_side),
                pow(1./3.,3)/3. * pow(1./4.,2)/2. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*y, OShape(0), front_side),
                pow(1./3.,3)/3. * pow(1./4.,2)/2. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(y*y*y*y*z, OShape(0), left_side),
                pow(1./4.,5)/5. * pow(1./5.,2)/2. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(y*y*y*y*z, OShape(0), right_side),
                pow(1./4.,5)/5. * pow(1./5.,2)/2. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*x*z, OShape(0), bottom_side),
                pow(1./3.,4)/4. * pow(1./5.,2)/2. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*x*z, OShape(0), top_side),
                pow(1./3.,4)/4. * pow(1./5.,2)/2. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*x*y*y*y*y, OShape(0), back_side),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*x*y*y*y*y, OShape(0), front_side),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * 1./6.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*x*y*y*y*y, OShape(0), early_side),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * 1./5.);
  assert_approx(rmesh3x4x5x6.intg_facerel_mon_on_oshape_side(x*x*x*y*y*y*y, OShape(0), late_side),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * 1./5.);
}

#[test]
fn test_intg_facerel_mon_x_facerel_poly_on_oshape_int() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_int(y, &poly(~[(2.,z*z), (3.,x*x*x*y*y*y*y*z)]), OShape(0)),
                2.*(1./3. * pow(1./4.,2)/2. * pow(1./5.,3)/3.) + 3.*(pow(1./3.,4)/4. * pow(1./4.,6)/6. * pow(1./5.,2)/2.));

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_int(x*z, &poly(~[(2.,y)]), OShape(0)),
                pow(1./3.,2)/2. * pow(1./4.,2) * pow(1./5.,2)/2.);

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_int(z*z, &poly(~[(2.,y), (3.,x*x)]), OShape(0)),
                2. * 1./3. * pow(1./4.,2)/2. * pow(1./5.,3)/3. + 3. * pow(1./3.,3)/3. * 1./4. * pow(1./5.,3)/3.);
}


#[test]
fn test_intg_facerel_mon_x_facerel_poly_on_oshape_side() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(one, &poly(~[(1.2,y*z*z), (3.4,x*x*x*y*y*y*y)]), OShape(0), left_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(one, &poly(~[(1.2,y*z*z), (3.4,x*x*x*y*y*y*y)]), OShape(0), right_side),
                1.2 * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(one, &poly(~[(1.5,y*z*z), (0.123,y*y*y*y*z)]), OShape(0), right_side),
                1.5 * pow(1./4.,2)/2. * pow(1./5.,3)/3. + 0.123 * pow(1./4.,5)/5. * pow(1./5.,2)/2.);

  // yz (2y + 4z) = 2y^2z + 4yz^2
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(y*z, &poly(~[(2.,y), (4.,z)]), OShape(0), right_side),
                2. * pow(1./4.,3)/3. * pow(1./5.,2)/2. 
              + 4. * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
 
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(y*z, &poly(~[(2.,y), (4.,z)]), OShape(0), left_side),
                2. * pow(1./4.,3)/3. * pow(1./5.,2)/2. 
              + 4. * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  
  // xz (2x + 4z) = 2x^2z + 4xz^2
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(x*z, &poly(~[(2.,x), (4.,z)]), OShape(0), top_side),
                2. * pow(1./3.,3)/3. * pow(1./5.,2)/2. 
              + 4. * pow(1./3.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(x*z, &poly(~[(2.,x), (4.,z)]), OShape(0), bottom_side),
                2. * pow(1./3.,3)/3. * pow(1./5.,2)/2. 
              + 4. * pow(1./3.,2)/2. * pow(1./5.,3)/3.);
  
  // xy (2x + 4y) = 2x^2y + 4xy^2
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(x*y, &poly(~[(2.,x), (4.,y)]), OShape(0), front_side),
                2. * pow(1./3.,3)/3. * pow(1./4.,2)/2. 
              + 4. * pow(1./3.,2)/2. * pow(1./4.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_poly_on_oshape_side(x*y, &poly(~[(2.,x), (4.,y)]), OShape(0), back_side),
                2. * pow(1./3.,3)/3. * pow(1./4.,2)/2. 
              + 4. * pow(1./3.,2)/2. * pow(1./4.,3)/3.);
}

#[test]
fn test_intg_intrel_mon_x_siderel_mon_on_oshape_side() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  // xy yz^2 = x y^2z^2 where x is 0 or 1/3 for the left and right sides, respectively.
  assert_approx(rmesh3x4x5.intg_intrel_mon_x_siderel_mon_on_oshape_side(x*y, y*z*z, OShape(0), left_side),
                0.);
  assert_approx(rmesh3x4x5.intg_intrel_mon_x_siderel_mon_on_oshape_side(x*y, y*z*z*z, OShape(0), right_side),
                1./3. * pow(1./4.,3)/3. * pow(1./5.,4)/4.);
  
  // yz xz^2 =  y xz^3 where y is 0 or 1/4 for the bottom and top sides, respectively.
  assert_approx(rmesh3x4x5.intg_intrel_mon_x_siderel_mon_on_oshape_side(y*z, x*z*z, OShape(0), bottom_side),
                0.);
  assert_approx(rmesh3x4x5.intg_intrel_mon_x_siderel_mon_on_oshape_side(y*z, x*z*z, OShape(0), top_side),
                1./4. * pow(1./3.,2)/2. * pow(1./5.,4)/4.);
  
  // yz xy^2 =  z xy^3 where z is 0 or 1/5 for back and front sides, respectively.
  assert_approx(rmesh3x4x5.intg_intrel_mon_x_siderel_mon_on_oshape_side(y*z, x*y*y, OShape(0), back_side),
                0.);
  assert_approx(rmesh3x4x5.intg_intrel_mon_x_siderel_mon_on_oshape_side(y*z, x*y*y, OShape(0), front_side),
                1./5. * pow(1./3.,2)/2. * pow(1./4.,4)/4.);
}

#[test]
fn test_intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side_dim0() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), left_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), right_side),
             pow(1./3.,3) * pow(1./4.,5)/5. * pow(1./5.,2)/2.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               y*z, &VectorMonomial{mon: x*x*x*y*y*y, mon_dim: Dim(0)}, OShape(0), right_side),
             pow(1./3.,3) * pow(1./4.,5)/5. * pow(1./5.,2)/2.);

  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               y*z, &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), left_side),
             -pow(1./4.,5)/5. * pow(1./5.,3)/3.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               y*z, &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(1)}, OShape(0), left_side),
             0.);
 
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               y*z, &VectorMonomial{mon: x*x*x*y, mon_dim: Dim(0)}, OShape(0), right_side),
             pow(1./3.,3) * pow(1./4.,3)/3. * pow(1./5.,2)/2.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               x*y*y*z, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), right_side),
             0.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), top_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), bottom_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), front_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), back_side),
             0.);
}  

fn test_intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side_dim1() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), back_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), front_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), left_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), right_side),
             0.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), bottom_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               x*z, &VectorMonomial{mon: x*x*x*z, mon_dim: Dim(1)}, OShape(0), bottom_side),
             pow(1./3.,5)/5. * pow(1./5.,3)/3.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), top_side),
             pow(1./4.,4) * pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               x*x*x*z, &VectorMonomial{mon: y*y*y*y, mon_dim: Dim(1)}, OShape(0), top_side),
             pow(1./4.,4) * pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               y*x*z, &VectorMonomial{mon: y*y*y*y, mon_dim: Dim(1)}, OShape(0), top_side),
             0.);
}

fn test_intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side_dim2() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), left_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), right_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), bottom_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), top_side),
             0.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               one, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), back_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               x*y, &VectorMonomial{mon: x*x*y*y*y*z, mon_dim: Dim(2)}, OShape(0), front_side),
             1./5. * pow(1./3.,4)/4. * pow(1./4.,5)/5.);

  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               x*y, &VectorMonomial{mon: x*x*y*y*y, mon_dim: Dim(2)}, OShape(0), back_side),
             pow(1./3.,4)/4. * pow(1./4.,5)/5.);
  assert_eq!(rmesh3x4x5.intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(
               x*y, &VectorMonomial{mon: x*x*y*y*y, mon_dim: Dim(2)}, OShape(0), front_side),
             pow(1./3.,4)/4. * pow(1./4.,5)/5.);
}


#[test]
fn test_intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side_dim0() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_poly = ~poly(~[(1.,one)]);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), left_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(0)}, OShape(0), right_side),
             pow(1./3.,3) * pow(1./4.,5)/5. * pow(1./5.,2)/2.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,y*z), (3.,y*y*z*z)]), &VectorMonomial{mon: x*x*x*y*y*y, mon_dim: Dim(0)}, OShape(0), right_side),
             2. * pow(1./3.,3) * pow(1./4.,5)/5. * pow(1./5.,2)/2.
           + 3. * pow(1./3.,3) * pow(1./4.,6)/6. * pow(1./5.,3)/3.);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,y*z), (3.,y*y*z*z)]), &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), left_side),
             -2. * pow(1./4.,5)/5. * pow(1./5.,3)/3.
           + -3. * pow(1./4.,6)/6. * pow(1./5.,4)/4.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y*z), (3.,x*y*y*z*z)]), &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), left_side),
             0.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y*z), (3.,x*y*y*z*z)]), &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), top_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y*z), (3.,x*y*y*z*z)]), &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), bottom_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y*z), (3.,x*y*y*z*z)]), &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), front_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y*z), (3.,x*y*y*z*z)]), &VectorMonomial{mon: y*y*y*z, mon_dim: Dim(0)}, OShape(0), back_side),
             0.);
}  

fn test_intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side_dim1() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_poly = ~poly(~[(1.,one)]);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), back_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), front_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), left_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), right_side),
             0.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), bottom_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*z), (3.,x*x*z*z)]), &VectorMonomial{mon: x*x*x*z, mon_dim: Dim(1)}, OShape(0), bottom_side),
             2. * pow(1./3.,5)/5. * pow(1./5.,3)/3.
           + 3. * pow(1./3.,6)/6. * pow(1./5.,4)/4.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(1)}, OShape(0), top_side),
             pow(1./4.,4) * pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*x*x*z), (-3.25,x*z*z)]), &VectorMonomial{mon: y*y*y*y, mon_dim: Dim(1)}, OShape(0), top_side),
             2. * pow(1./4.,4) * pow(1./3.,4)/4. * pow(1./5.,2)/2.
           - 3.25 * pow(1./4.,4) * pow(1./3.,3)/3. * pow(1./5.,3)/3.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,y*x*z), (5.2,y*z*z)]), &VectorMonomial{mon: y*y*y*y, mon_dim: Dim(1)}, OShape(0), top_side),
             0.);
}

fn test_intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side_dim2() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1f64, 2., 3.],
                                                   ~[2f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));
  let back_side = lesser_side_face_perp_to_axis(Dim(2));
  let front_side = greater_side_face_perp_to_axis(Dim(2));

  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_poly = ~poly(~[(1.,one)]);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), left_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), right_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), bottom_side),
             0.);
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), top_side),
             0.);
  
  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               one_poly, &VectorMonomial{mon: x*x*x*y*y*y*y*z, mon_dim: Dim(2)}, OShape(0), back_side),
             0.);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y), (5.2,x*x*y*y)]), &VectorMonomial{mon: x*x*y*y*y*z, mon_dim: Dim(2)}, OShape(0), front_side),
             2. * 1./5. * pow(1./3.,4)/4. * pow(1./4.,5)/5.
          + 5.2 * 1./5. * pow(1./3.,5)/5. * pow(1./4.,6)/6.);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y), (5.2,x*x*y*y)]), &VectorMonomial{mon: x*x*y*y*y, mon_dim: Dim(2)}, OShape(0), back_side),
             2. * pow(1./3.,4)/4. * pow(1./4.,5)/5.
         + 5.2  * pow(1./3.,5)/5. * pow(1./4.,6)/6.);

  assert_eq!(rmesh3x4x5.intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side(
               &poly(~[(2.,x*y), (5.2,x*x*y*y)]), &VectorMonomial{mon: x*x*y*y*y, mon_dim: Dim(2)}, OShape(0), front_side),
             2. * pow(1./3.,4)/4. * pow(1./4.,5)/5.
         + 5.2  * pow(1./3.,5)/5. * pow(1./4.,6)/6.);
}


fn assert_approx(a:R, b:R) -> () {
  assert!(abs(a - b) < 10e-9)
}


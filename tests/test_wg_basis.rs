use wg_basis::{WGBasis, BasisElNum, FaceMonNum};
use common::{Deg, Dim};
use mesh::{FENum, OShape, SideFace, NBSideNum, NBSideInclusions};
use rectangle_mesh::{RectMesh, MeshCoord};
use monomial::{Mon2d, MaxMonDeg};


/*
 3 cols x 2 rows mesh, k = 2
 ----------
 |  |  |  |
 ----------
 |  |  |  |
 ----------

 6 interior monomials, 2 side monomials
 interior basis els: 3 * 2 * 6 = 36
 vertical side basis els: 2 * 2 * 2 = 8
 horizontal side basis els: 3 * 2 = 6
*/

#[test]
fn test_int_mons_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(basis.ref_int_mons(), [one, y, y*y, x, x*y, x*x]);

  assert_eq!(basis.mons_per_fe_int(), 6);
}


#[test]
fn test_side_mons_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(basis.side_mons_by_dep_dim[0], ~[one, y]);
  assert_eq!(basis.side_mons_by_dep_dim[1], ~[one, x]);
  
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(0)), [one, y]); // left side
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(1)), [one, y]); // right side
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(2)), [one, x]); // bottom side
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(3)), [one, x]); // top side
  
  assert_eq!(basis.mons_per_fe_side(), 2);
}


#[test]
fn test_is_int_supp_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert!(basis.is_int_supported(BasisElNum(0)));
  assert!(basis.is_int_supported(BasisElNum(35)));
  assert!(!basis.is_int_supported(BasisElNum(36)));
}

#[test]
fn test_first_nb_side_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.first_nb_side_beln, BasisElNum(36));
}

#[test]
fn test_nb_side_num_blocks_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));

  assert_eq!(basis.mesh.perp_axis_for_nb_side(basis.support_nb_side_num(BasisElNum(36))), Dim(0));
  assert_eq!(basis.mesh.perp_axis_for_nb_side(basis.support_nb_side_num(BasisElNum(43))), Dim(0));
  assert_eq!(basis.mesh.perp_axis_for_nb_side(basis.support_nb_side_num(BasisElNum(44))), Dim(1));
  assert_eq!(basis.mesh.perp_axis_for_nb_side(basis.support_nb_side_num(BasisElNum(49))), Dim(1));
}

#[test]
#[should_fail]
fn test_bad_nb_side_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  basis.mesh.perp_axis_for_nb_side(basis.support_nb_side_num(BasisElNum(50)));
}

// The first section of basis elements are monomials which are assigned to interiors in blocks of 6.
#[test]
fn test_int_support_fes_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));

  assert_eq!(basis.support_int_fe_num(BasisElNum(0)), FENum(0));
  assert_eq!(basis.support_int_fe_num(BasisElNum(5)), FENum(0));
  assert_eq!(basis.support_int_fe_num(BasisElNum(6)), FENum(1));
  assert_eq!(basis.support_int_fe_num(BasisElNum(11)), FENum(1));
  assert_eq!(basis.support_int_fe_num(BasisElNum(12)), FENum(2));
  assert_eq!(basis.support_int_fe_num(BasisElNum(29)), FENum(4));
  assert_eq!(basis.support_int_fe_num(BasisElNum(30)), FENum(5));
  assert_eq!(basis.support_int_fe_num(BasisElNum(35)), FENum(5));
}

// The next section of basis elements consists of those non-boundary sides perpendicular to the x axis.
#[test]
fn test_vert_nb_side_fe_incls_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));

  let left_face = SideFace(0);
  let right_face = SideFace(1);

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(36)),
             NBSideInclusions { nb_side_num: NBSideNum(0),
                                fe1: FENum(0),
                                side_face_in_fe1: right_face,
                                fe2: FENum(1),
                                side_face_in_fe2: left_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(37)),
             NBSideInclusions { nb_side_num: NBSideNum(0),
                                fe1: FENum(0),
                                side_face_in_fe1: right_face,
                                fe2: FENum(1),
                                side_face_in_fe2: left_face });
  
  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(38)),
             NBSideInclusions { nb_side_num: NBSideNum(1),
                                fe1: FENum(1),
                                side_face_in_fe1: right_face,
                                fe2: FENum(2),
                                side_face_in_fe2: left_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(39)),
             NBSideInclusions { nb_side_num: NBSideNum(1),
                                fe1: FENum(1),
                                side_face_in_fe1: right_face,
                                fe2: FENum(2),
                                side_face_in_fe2: left_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(40)),
             NBSideInclusions { nb_side_num: NBSideNum(2),
                                fe1: FENum(3),
                                side_face_in_fe1: right_face,
                                fe2: FENum(4),
                                side_face_in_fe2: left_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(41)),
             NBSideInclusions { nb_side_num: NBSideNum(2),
                                fe1: FENum(3),
                                side_face_in_fe1: right_face,
                                fe2: FENum(4),
                                side_face_in_fe2: left_face });
  
  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(42)),
             NBSideInclusions { nb_side_num: NBSideNum(3),
                                fe1: FENum(4),
                                side_face_in_fe1: right_face,
                                fe2: FENum(5),
                                side_face_in_fe2: left_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(43)),
             NBSideInclusions { nb_side_num: NBSideNum(3),
                                fe1: FENum(4),
                                side_face_in_fe1: right_face,
                                fe2: FENum(5),
                                side_face_in_fe2: left_face });
}

// The final section of basis elements contains the horizontal side monomials. 
#[test]
fn test_horiz_nb_side_fe_incls_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));

  let bottom_face = SideFace(2);
  let top_face = SideFace(3);

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(44)),
             NBSideInclusions { nb_side_num: NBSideNum(4),
                                fe1: FENum(0),
                                side_face_in_fe1: top_face,
                                fe2: FENum(3),
                                side_face_in_fe2: bottom_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(45)),
             NBSideInclusions { nb_side_num: NBSideNum(4),
                                fe1: FENum(0),
                                side_face_in_fe1: top_face,
                                fe2: FENum(3),
                                side_face_in_fe2: bottom_face });
  
  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(46)),
             NBSideInclusions { nb_side_num: NBSideNum(5),
                                fe1: FENum(1),
                                side_face_in_fe1: top_face,
                                fe2: FENum(4),
                                side_face_in_fe2: bottom_face });

  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(47)),
             NBSideInclusions { nb_side_num: NBSideNum(5),
                                fe1: FENum(1),
                                side_face_in_fe1: top_face,
                                fe2: FENum(4),
                                side_face_in_fe2: bottom_face });
  
  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(48)),
             NBSideInclusions { nb_side_num: NBSideNum(6),
                                fe1: FENum(2),
                                side_face_in_fe1: top_face,
                                fe2: FENum(5),
                                side_face_in_fe2: bottom_face });
  
  assert_eq!(basis.fe_inclusions_of_side_support(BasisElNum(49)),
             NBSideInclusions { nb_side_num: NBSideNum(6),
                                fe1: FENum(2),
                                side_face_in_fe1: top_face,
                                fe2: FENum(5),
                                side_face_in_fe2: bottom_face });
}

// interior monomials, first finite element

#[test]
fn test_fe0_int_mon_retrieval_by_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(basis.int_mon(BasisElNum(0)), one);
  assert_eq!(basis.int_mon(BasisElNum(1)), y);
  assert_eq!(basis.int_mon(BasisElNum(2)), y*y);
  assert_eq!(basis.int_mon(BasisElNum(3)), x);
  assert_eq!(basis.int_mon(BasisElNum(4)), x*y);
  assert_eq!(basis.int_mon(BasisElNum(5)), x*x);
}

#[test]
fn test_fe0_int_rel_mon_num_retrieval_by_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.int_rel_mon_num(BasisElNum(0)), FaceMonNum(0));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(1)), FaceMonNum(1));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(2)), FaceMonNum(2));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(3)), FaceMonNum(3));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(4)), FaceMonNum(4));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(5)), FaceMonNum(5));
}


#[test]
fn test_fe0_int_mon_beln_by_fe_and_face_mon_num_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.int_mon_el_num(FENum(0), FaceMonNum(0)), BasisElNum(0));
  assert_eq!(basis.int_mon_el_num(FENum(0), FaceMonNum(1)), BasisElNum(1));
  assert_eq!(basis.int_mon_el_num(FENum(0), FaceMonNum(2)), BasisElNum(2));
  assert_eq!(basis.int_mon_el_num(FENum(0), FaceMonNum(3)), BasisElNum(3));
  assert_eq!(basis.int_mon_el_num(FENum(0), FaceMonNum(4)), BasisElNum(4));
  assert_eq!(basis.int_mon_el_num(FENum(0), FaceMonNum(5)), BasisElNum(5));
}

// interior monomials, second finite element

#[test]
fn test_fe1_int_mon_retrieval_by_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(basis.int_mon(BasisElNum(6)), one);
  assert_eq!(basis.int_mon(BasisElNum(7)), y);
  assert_eq!(basis.int_mon(BasisElNum(8)), y*y);
  assert_eq!(basis.int_mon(BasisElNum(9)), x);
  assert_eq!(basis.int_mon(BasisElNum(10)), x*y);
  assert_eq!(basis.int_mon(BasisElNum(11)), x*x);
}

#[test]
fn test_fe1_int_rel_mon_num_retrieval_by_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.int_rel_mon_num(BasisElNum(6)), FaceMonNum(0));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(7)), FaceMonNum(1));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(8)), FaceMonNum(2));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(9)), FaceMonNum(3));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(10)), FaceMonNum(4));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(11)), FaceMonNum(5));
}


#[test]
fn test_fe1_int_mon_beln_by_fe_and_face_mon_num_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.int_mon_el_num(FENum(1), FaceMonNum(0)), BasisElNum(6));
  assert_eq!(basis.int_mon_el_num(FENum(1), FaceMonNum(1)), BasisElNum(7));
  assert_eq!(basis.int_mon_el_num(FENum(1), FaceMonNum(2)), BasisElNum(8));
  assert_eq!(basis.int_mon_el_num(FENum(1), FaceMonNum(3)), BasisElNum(9));
  assert_eq!(basis.int_mon_el_num(FENum(1), FaceMonNum(4)), BasisElNum(10));
  assert_eq!(basis.int_mon_el_num(FENum(1), FaceMonNum(5)), BasisElNum(11));
}

// interior monomials, last finite element

#[test]
fn test_fe5_int_mon_retrieval_by_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(basis.int_mon(BasisElNum(30)), one);
  assert_eq!(basis.int_mon(BasisElNum(31)), y);
  assert_eq!(basis.int_mon(BasisElNum(32)), y*y);
  assert_eq!(basis.int_mon(BasisElNum(33)), x);
  assert_eq!(basis.int_mon(BasisElNum(34)), x*y);
  assert_eq!(basis.int_mon(BasisElNum(35)), x*x);
}

#[test]
fn test_fe5_int_rel_mon_num_retrieval_by_beln_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.int_rel_mon_num(BasisElNum(30)), FaceMonNum(0));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(31)), FaceMonNum(1));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(32)), FaceMonNum(2));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(33)), FaceMonNum(3));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(34)), FaceMonNum(4));
  assert_eq!(basis.int_rel_mon_num(BasisElNum(35)), FaceMonNum(5));
}


#[test]
fn test_fe5_int_mon_beln_by_fe_and_face_mon_num_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.int_mon_el_num(FENum(5), FaceMonNum(0)), BasisElNum(30));
  assert_eq!(basis.int_mon_el_num(FENum(5), FaceMonNum(1)), BasisElNum(31));
  assert_eq!(basis.int_mon_el_num(FENum(5), FaceMonNum(2)), BasisElNum(32));
  assert_eq!(basis.int_mon_el_num(FENum(5), FaceMonNum(3)), BasisElNum(33));
  assert_eq!(basis.int_mon_el_num(FENum(5), FaceMonNum(4)), BasisElNum(34));
  assert_eq!(basis.int_mon_el_num(FENum(5), FaceMonNum(5)), BasisElNum(35));
}

// vertical side monomials

#[test]
fn test_side_mons_seq_between_fe0_and_fe1_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let left_face = SideFace(0);
  let right_face = SideFace(1);
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  let fe0_right_side_mons = basis.side_mons_for_fe_side(FENum(0), right_face);
  assert_eq!(fe0_right_side_mons, [one, y]);
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), right_face), fe0_right_side_mons);
  assert_eq!(fe0_right_side_mons, basis.side_mons_for_fe_side(FENum(1), left_face));
  assert_eq!(fe0_right_side_mons, basis.side_mons_for_oshape_side(OShape(0), left_face));
}

#[test]
fn test_side_mon_beln_between_fe0_and_fe1_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let left_face = SideFace(0);
  let right_face = SideFace(1);
  
  assert_eq!(basis.fe_side_mon_el_num(FENum(0), right_face, FaceMonNum(0)), BasisElNum(36));
  assert_eq!(basis.fe_side_mon_el_num(FENum(1), left_face,  FaceMonNum(0)), BasisElNum(36));

  assert_eq!(basis.fe_side_mon_el_num(FENum(0), right_face, FaceMonNum(1)), BasisElNum(37));
  assert_eq!(basis.fe_side_mon_el_num(FENum(1), left_face,  FaceMonNum(1)), BasisElNum(37));
}

#[test]
fn test_side_rel_mon_nums_by_beln_between_fe0_and_fe1_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.side_rel_mon_num(BasisElNum(36)), FaceMonNum(0));
  assert_eq!(basis.side_rel_mon_num(BasisElNum(37)), FaceMonNum(1));
}

#[test]
fn test_side_mons_seq_between_fe1_and_fe2_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let left_face = SideFace(0);
  let right_face = SideFace(1);
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  let fe1_right_side_mons = basis.side_mons_for_fe_side(FENum(1), right_face);
  assert_eq!(fe1_right_side_mons, [one, y]);
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), right_face), fe1_right_side_mons);
  assert_eq!(fe1_right_side_mons, basis.side_mons_for_fe_side(FENum(1), left_face));
  assert_eq!(fe1_right_side_mons, basis.side_mons_for_oshape_side(OShape(0), left_face));
}

#[test]
fn test_side_mon_beln_between_fe1_and_fe2_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let left_face = SideFace(0);
  let right_face = SideFace(1);
  
  assert_eq!(basis.fe_side_mon_el_num(FENum(1), right_face, FaceMonNum(0)), BasisElNum(38));
  assert_eq!(basis.fe_side_mon_el_num(FENum(2), left_face,  FaceMonNum(0)), BasisElNum(38));

  assert_eq!(basis.fe_side_mon_el_num(FENum(1), right_face, FaceMonNum(1)), BasisElNum(39));
  assert_eq!(basis.fe_side_mon_el_num(FENum(2), left_face,  FaceMonNum(1)), BasisElNum(39));
}

#[test]
fn test_side_rel_mon_nums_by_beln_between_fe1_and_fe2_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.side_rel_mon_num(BasisElNum(38)), FaceMonNum(0));
  assert_eq!(basis.side_rel_mon_num(BasisElNum(39)), FaceMonNum(1));
}

#[test]
fn test_side_mon_beln_between_fe4_and_fe5_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let left_face = SideFace(0);
  let right_face = SideFace(1);

  assert_eq!(basis.fe_side_mon_el_num(FENum(4), right_face, FaceMonNum(0)), BasisElNum(42));
  assert_eq!(basis.fe_side_mon_el_num(FENum(5), left_face,  FaceMonNum(0)), BasisElNum(42));

  assert_eq!(basis.fe_side_mon_el_num(FENum(4), right_face, FaceMonNum(1)), BasisElNum(43));
  assert_eq!(basis.fe_side_mon_el_num(FENum(5), left_face,  FaceMonNum(1)), BasisElNum(43));
}

#[test]
fn test_side_rel_mon_nums_by_beln_between_fe4_and_fe5_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.side_rel_mon_num(BasisElNum(42)), FaceMonNum(0));
  assert_eq!(basis.side_rel_mon_num(BasisElNum(43)), FaceMonNum(1));
}


// horizontal side monomials

#[test]
fn test_side_mons_seq_between_fe0_and_fe3_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let bottom_face = SideFace(2);
  let top_face = SideFace(3);
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };

  let fe0_top_side_mons = basis.side_mons_for_fe_side(FENum(0), top_face);
  assert_eq!(fe0_top_side_mons, [one, x]);
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), top_face), fe0_top_side_mons);
  assert_eq!(fe0_top_side_mons, basis.side_mons_for_fe_side(FENum(3), bottom_face));
  assert_eq!(fe0_top_side_mons, basis.side_mons_for_oshape_side(OShape(0), bottom_face));
}

#[test]
fn test_side_mon_beln_between_fe0_and_fe3_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let bottom_face = SideFace(2);
  let top_face = SideFace(3);
  
  assert_eq!(basis.fe_side_mon_el_num(FENum(0), top_face,    FaceMonNum(0)), BasisElNum(44));
  assert_eq!(basis.fe_side_mon_el_num(FENum(3), bottom_face, FaceMonNum(0)), BasisElNum(44));

  assert_eq!(basis.fe_side_mon_el_num(FENum(0), top_face,    FaceMonNum(1)), BasisElNum(45));
  assert_eq!(basis.fe_side_mon_el_num(FENum(3), bottom_face, FaceMonNum(1)), BasisElNum(45));
}

#[test]
fn test_side_rel_mon_nums_by_beln_between_fe0_and_fe3_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  assert_eq!(basis.side_rel_mon_num(BasisElNum(44)), FaceMonNum(0));
  assert_eq!(basis.side_rel_mon_num(BasisElNum(45)), FaceMonNum(1));
}

#[test]
fn test_wgrads_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));

  let xy_on_int_wgrad = basis.int_mon_wgrad(FaceMonNum(4), OShape(0));
  assert_eq!(xy_on_int_wgrad.comp_mon_coefs[0].as_slice(), &[3./2., 0., -3.]); // See test_weak_gradient.rs.
  assert_eq!(xy_on_int_wgrad.comp_mon_coefs[1].as_slice(), &[3./2., -3., 0.]); // 

  let y_on_right_side_wgrad = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(1));
  assert_eq!(y_on_right_side_wgrad.comp_mon_coefs[0].as_slice(), &[-3./2., 1., 3.]);
  assert_eq!(y_on_right_side_wgrad.comp_mon_coefs[1].as_slice(), &[0., 0., 0.]);

  let x_on_top_side_wgrad = basis.side_mon_wgrad(FaceMonNum(1), OShape(0), SideFace(3));
  assert_eq!(x_on_top_side_wgrad.comp_mon_coefs[0].as_slice(), &[0., 0., 0.]);
  assert_eq!(x_on_top_side_wgrad.comp_mon_coefs[1].as_slice(), &[-3./2., 3., 1.]);
}

#[test]
fn test_int_mons_3x2_deg3() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  
  let basis = WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  
  assert_eq!(basis.ref_int_mons(), [one, y, y*y, y*y*y, x, x*y, x*y*y, x*x, x*x*y, x*x*x]);
  
  assert_eq!(basis.mons_per_fe_int(), 10);
}

#[test]
fn test_side_mons_3x2_deg3() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  
  let basis = WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));
  
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(basis.side_mons_by_dep_dim[0], ~[one, y, y*y]);
  assert_eq!(basis.side_mons_by_dep_dim[1], ~[one, x, x*x]);
  
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(0)), [one, y, y*y]); // left side
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(1)), [one, y, y*y]); // right side
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(2)), [one, x, x*x]); // bottom side
  assert_eq!(basis.side_mons_for_oshape_side(OShape(0), SideFace(3)), [one, x, x*x]); // top side

  assert_eq!(basis.mons_per_fe_side(), 3);
}


#[test]
fn test_interacting_els_est_3x2_deg3() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(3), MaxMonDeg(2));

  let mut int_int_inters = 0u;
  let mut int_side_and_vv_inters = 0u;
  let mut side_side_inters = 0u;

  for el_1 in range(0, basis.num_els()) {
    for el_2 in range(0, basis.num_els()) {
      let (el_1, el_2) = (BasisElNum(el_1), BasisElNum(el_2));

      // both elements interior supported
      if basis.is_int_supported(el_1) && basis.is_int_supported(el_2) {
        let (fe_1, fe_2) = (basis.support_int_fe_num(el_1), basis.support_int_fe_num(el_2));
        if fe_1 == fe_2 {
          int_int_inters += 1;
        }
      }

      // interior - side
      else if basis.is_int_supported(el_1) && basis.is_side_supported(el_2) {
        let (fe_1, incls_2) = (basis.support_int_fe_num(el_1), basis.fe_inclusions_of_side_support(el_2));
        if fe_1 == incls_2.fe1 || fe_1 == incls_2.fe2 {
          int_side_and_vv_inters += 1;
        }
      }
      
      // side - interior
      else if basis.is_side_supported(el_1) && basis.is_int_supported(el_2) {
        let (incls_1, fe_2) = (basis.fe_inclusions_of_side_support(el_1), basis.support_int_fe_num(el_2));
        if fe_2 == incls_1.fe1 || fe_2 == incls_1.fe2 {
          int_side_and_vv_inters += 1;
        }
      }
      
      // side - side 
      else if basis.is_side_supported(el_1) && basis.is_side_supported(el_2) {
        let (incls_1, incls_2) = (basis.fe_inclusions_of_side_support(el_1), basis.fe_inclusions_of_side_support(el_2));

        if incls_1.fe1 == incls_2.fe1 || incls_1.fe2 == incls_2.fe2 || incls_1.fe1 == incls_2.fe2 || incls_1.fe2 == incls_2.fe1 {
          side_side_inters += 1;
        }
      }
      else { fail!("Support for basis elements did not match exhaustive alternatives."); }
    }
  }
  
  assert_eq!(int_int_inters + int_side_and_vv_inters + side_side_inters,
             basis.est_num_el_el_pairs_with_common_supp_fes(false));
}

#[test]
fn test_interacting_els_est_5x6_deg4() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(5),MeshCoord(6)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(4), MaxMonDeg(3));

  let mut int_int_inters = 0u;
  let mut int_side_and_vv_inters = 0u;
  let mut side_side_inters = 0u;
  let mut upper_triangle_inters = 0u;

  for el_1 in range(0, basis.num_els()) {
    for el_2 in range(0, basis.num_els()) {
      let (el_1, el_2) = (BasisElNum(el_1), BasisElNum(el_2));

      // both elements interior supported
      if basis.is_int_supported(el_1) && basis.is_int_supported(el_2) {
        let (fe_1, fe_2) = (basis.support_int_fe_num(el_1), basis.support_int_fe_num(el_2));
        if fe_1 == fe_2 {
          int_int_inters += 1;
          if *el_1 <= *el_2 { upper_triangle_inters += 1; }
        }
      }

      // interior - side
      else if basis.is_int_supported(el_1) && basis.is_side_supported(el_2) {
        let (fe_1, incls_2) = (basis.support_int_fe_num(el_1), basis.fe_inclusions_of_side_support(el_2));
        if fe_1 == incls_2.fe1 || fe_1 == incls_2.fe2 {
          int_side_and_vv_inters += 1;
          if *el_1 <= *el_2 { upper_triangle_inters += 1; }
        }
      }
      
      // side - interior
      else if basis.is_side_supported(el_1) && basis.is_int_supported(el_2) {
        let (incls_1, fe_2) = (basis.fe_inclusions_of_side_support(el_1), basis.support_int_fe_num(el_2));
        if fe_2 == incls_1.fe1 || fe_2 == incls_1.fe2 {
          int_side_and_vv_inters += 1;
          if *el_1 <= *el_2 { upper_triangle_inters += 1; }
        }
      }
      
      // side - side 
      else if basis.is_side_supported(el_1) && basis.is_side_supported(el_2) {
        let (incls_1, incls_2) = (basis.fe_inclusions_of_side_support(el_1), basis.fe_inclusions_of_side_support(el_2));

        if incls_1.fe1 == incls_2.fe1 || incls_1.fe2 == incls_2.fe2 || incls_1.fe1 == incls_2.fe2 || incls_1.fe2 == incls_2.fe1 {
          side_side_inters += 1;
          if *el_1 <= *el_2 { upper_triangle_inters += 1; }
        }
      }
      else { fail!("Support for basis elements did not match exhaustive alternatives."); }
    }
  }
 
  let tot_inters = int_int_inters + int_side_and_vv_inters + side_side_inters;

  assert_eq!(tot_inters,
             basis.est_num_el_el_pairs_with_common_supp_fes(false));

  assert_eq!(basis.est_num_el_el_pairs_with_common_supp_fes(true), upper_triangle_inters);
}

#[test]
fn test_int_L2_inner_products_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let int_ips = basis.ips_int_mons_for_oshape(OShape(0));
  assert_eq!(int_ips.get(0,0), 1.);    // one vs one
  assert_eq!(int_ips.get(0,1), 1./2.); // one vs y
  assert_eq!(int_ips.get(0,2), 1./3.); // one vs y^2
  assert_eq!(int_ips.get(0,3), 1./2.); // one vs x
  assert_eq!(int_ips.get(0,4), 1./4.); // one vs xy
  assert_eq!(int_ips.get(0,5), 1./3.); // one vs x^2
}

#[test]
fn test_side_L2_inner_products_3x2_deg2() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[3.,2.], ~[MeshCoord(3),MeshCoord(2)]);
  let basis = WGBasis::new(rmesh, MaxMonDeg(2), MaxMonDeg(1));
  
  let right_ips = basis.ips_side_mons_for_oshape_side(OShape(0), SideFace(3));
  assert_eq!(right_ips.get(0,0), 1.);    // one vs one
  assert_eq!(right_ips.get(0,1), 1./2.); // one vs y
  assert_eq!(right_ips.get(1,1), 1./3.); // y vs y

  let top_ips = basis.ips_side_mons_for_oshape_side(OShape(0), SideFace(1));
  assert_eq!(top_ips.get(0,0), 1.);    // one vs one
  assert_eq!(top_ips.get(0,1), 1./2.); // one vs x
  assert_eq!(top_ips.get(1,1), 1./3.); // x vs x
}


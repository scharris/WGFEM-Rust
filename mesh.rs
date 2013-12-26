use common::*;
use polynomial::Polynomial;
use vector_monomial::VectorMonomial;

// An FENum identifies a finite element in the mesh.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone,IterBytes)]
pub struct FENum(uint);

// An NBSideNum identifies a mesh side among all sides in the mesh.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone,IterBytes)]
pub struct NBSideNum(uint);

// A SideFace identifies a side within the context of a single oriented shape
// or finite element.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone,IterBytes)]
pub struct SideFace(uint);

/// An identifier for the shape of a finite element together with the orientation
/// of the shape (rotation will yield a different oriented shape). Many calculations
/// on a finite element can be expressed in element-local coordinates and will only
/// depend on the oriented shape of the element.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct OShape(uint);


/// The representations of a non-boundary side as a pair of finite element faces.
/// The elements are arranged so that fe1 < fe2.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct NBSideInclusions {
  nb_side_num: NBSideNum,

  fe1: FENum,
  side_face_in_fe1: SideFace,

  fe2: FENum,
  side_face_in_fe2: SideFace
}
impl NBSideInclusions {
  pub fn new(nb_side_num: NBSideNum, fe1: FENum, side_face_in_fe1: SideFace, fe2: FENum, side_face_in_fe2: SideFace) -> NBSideInclusions
  {
    NBSideInclusions{ nb_side_num: nb_side_num,
                      fe1: fe1, side_face_in_fe1: side_face_in_fe1,
                      fe2: fe2, side_face_in_fe2: side_face_in_fe2 }
  }
}

pub trait Mesh<Mon> {

  fn num_fes(&self) -> uint;

  fn num_nb_sides(&self) -> uint;

  fn num_oriented_element_shapes(&self) -> uint;

  fn oriented_shape_for_fe(&self, fe: FENum) -> OShape;

  fn num_side_faces_for_oshape(&self, os: OShape) -> uint;

  // The dependent dim is, for a given side, a chosen dimension j which is affine-dependent on
  // the other dimensions on the side. That is, the function returns a j for which c_0,...,c_d
  // exist such that
  //   x_j = c_0 + sum_{i=1..d} c_i x_i for all (x_1,...,x_d) in the side.
  // There may be more than one such coordinate number, in which case any one of these may be
  // returned depending on the mesh implementation.
  fn dependent_dim_for_oshape_side(&self, os: OShape, sf: SideFace) -> Dim;

  fn fe_inclusions_of_nb_side(&self, side_num: NBSideNum) -> NBSideInclusions;

  // Return non-boundary side number of the indicated fe relative side.
  fn nb_side_num_for_fe_side(&self, fe: FENum, sf: SideFace) -> NBSideNum;

  fn is_boundary_side(&self, fe: FENum, sf: SideFace) -> bool;

  fn num_boundary_sides(&self) -> uint;
  
  fn boundary_fes_by_oshape_side(&self) -> ~[~[~[FENum]]]; // oshape, side face -> fes
  
  fn shape_diameter_inv(&self, os: OShape) -> R;

  fn max_fe_diameter(&self) -> R;

  fn num_nb_sides_for_fe(&self, fe: FENum) -> uint;

  fn max_num_shape_sides(&self) -> uint;

  // integration functions
  
  fn intg_global_fn_on_fe_int(&self, f: |&[R]| -> R, fe: FENum) -> R;

  fn intg_global_fn_x_facerel_mon_on_fe_int(&self, f: |&[R]| -> R, mon: Mon, fe: FENum) -> R;

  fn intg_global_fn_x_facerel_mon_on_fe_side(&self, f: |&[R]| -> R, mon: Mon, fe: FENum, sf: SideFace) -> R;
  
  fn intg_mixed_global_and_facerel_fn_on_fe_int(&self, f: |&[R], &[R]| -> R, fe: FENum) -> R; 

  fn intg_facerel_poly_on_oshape_int<P:Polynomial<Mon>>(&self, p: &P, os: OShape) -> R;

  fn intg_facerel_poly_x_facerel_poly_on_oshape_side<P:Polynomial<Mon>>(&self, p1: &P, p2: &P, os: OShape, sf: SideFace) -> R;

  fn intg_facerel_mon_on_oshape_int(&self, mon: Mon, os: OShape) -> R;

  fn intg_facerel_mon_on_oshape_side(&self, mon: Mon, os: OShape, sf: SideFace) -> R;

  fn intg_facerel_mon_x_facerel_poly_on_oshape_side<P:Polynomial<Mon>>(&self, mon: Mon, p: &P, os: OShape, sf: SideFace) -> R;

  fn intg_intrel_mon_x_siderel_mon_on_oshape_side(&self, int_mon: Mon, side_mon: Mon, os: OShape, sf: SideFace) -> R;
  
  fn intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(&self, mon: Mon, q: &VectorMonomial<Mon>, os: OShape, sf: SideFace) -> R;
 
}


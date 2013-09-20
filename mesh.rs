extern mod extra;
use extra::treemap::TreeMap;
use std::vec;
use std::iter::Zip;
use common::*;
use monomial::Monomial;
use polynomial::Polynomial;
use vector_monomial::VectorMonomial;

mod common;
mod monomial;
mod polynomial;
mod vector_monomial;


// An FENum identifies a finite element in the mesh.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct FENum(uint);

// An NBSideNum identifies a mesh side among all sides in the mesh.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct NBSideNum(uint);

// A SideFace identifies a side within the context of a single oriented shape
// or finite element.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct SideFace(u8);

// A Face locally identifies either the interior or a side in the
// context of a single oriented shape or finite element.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub enum Face {
  Interior,
  Side(SideFace)
}

// An identifier for the shape of a finite element together with the orientation
// of the shape (rotation will yield a different oriented shape). Many calculations
// on a finite element can be expressed in element-local coordinates and will only
// depend on the oriented shape of the element.
#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct OShape(u32);


#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct NBSideInclusions {
  nb_side_num: NBSideNum,

  fe1: FENum,
  sideface_in_fe1: SideFace,

  fe2: FENum,
  sideface_in_fe2: SideFace
}

pub trait Mesh<M> {

  fn num_fes(&self) -> uint;

  fn num_nb_sides(&self) -> uint;

  fn num_oriented_element_shapes(&self) -> uint;

  fn oriented_shape_for_fe(&self, fe: FENum) -> OShape;

  fn num_side_faces_for_fe(&self, fe: FENum) -> uint;

  fn num_side_faces_for_shape(&self, oshape: OShape) -> uint;

  // The dependent dim is, for a given side, the dimension j which is affine-
  // dependent on the other dimensions on the side. That is, the function
  // returns a j for which c_0,...,c_d exist such that
  //   x_j = c_0 + sum_{i=1..d} c_i x_i for all (x_1,...,x_d) in the side.
  // There may be more than one such coordinate number, in which case any 
  // one of these may be returned depending on the mesh implementation.
  fn dependent_dim_for_oshape_side(&self, oshape: OShape, side_face: SideFace) -> Dim;

  fn fe_inclusions_of_nb_side(&self, side_num: NBSideNum) -> NBSideInclusions;

  // Return non-boundary side number of the indicated fe relative side.
  fn nb_side_num_for_fe_side(&self, fe: FENum, side_face: SideFace) -> NBSideNum;

  fn is_boundary_side(&self, fe: FENum, side_face: SideFace) -> bool;

  fn num_boundary_sides(&self) -> uint;
  
  fn shape_diameter_inv(&self, oshape: OShape) -> R;

  fn max_fe_diameter(&self) -> R;

  fn fe_interior_origin(&self, fe: FENum) -> ~[R];

  fn num_non_boundary_sides_for_fe(&self, fe: FENum) -> uint;

  fn max_num_shape_sides(&self) -> uint;


  // integration functions
  
  fn intg_global_fn_on_fe_int(&self,
       g: &fn(&[R]) -> R,
       fe: FENum) -> R;
  
  fn intg_global_fn_x_facerel_mon_on_fe_int(&self,
       g: &fn(&[R]) -> R,
       mon: M,
       fe: FENum) -> R;
 
  fn intg_facerel_poly_on_oshape_int<P:Polynomial<M>>(&self,
       p: &P,
       oshape: OShape) -> R;

  fn intg_facerel_poly_x_facerel_poly_on_oshape_int<P:Polynomial<M>>(&self,
       p1: &P,
       p2: &P,
       oshape: OShape) -> R;

  fn intg_facerel_poly_x_facerel_poly_on_oshape_side<P:Polynomial<M>>(&self,
       p1: &P,
       p2: &P,
       oshape: OShape,
       side_face: SideFace) -> R;

  fn intg_facerel_mon_on_oshape_int(&self,
       mon: M,
       oshape: OShape) -> R;

  fn intg_facerel_mon_on_oshape_side(&self,
       mon: M,
       oshape: OShape,
       side_face: SideFace) -> R;

  fn intg_facerel_mon_x_facerel_poly_on_oshape_int<P:Polynomial<M>>(&self,
       mon: M,
       p: &P,
       oshape: OShape) -> R;

  fn intg_facerel_mon_x_facerel_poly_on_oshape_side<P:Polynomial<M>>(&self,
       mon: M,
       p: &P,
       oshape: OShape,
       side_face: SideFace) -> R;

  fn intg_intrel_mon_x_siderel_mon_on_oshape_side(&self,
       int_mon: M,
       side_mon: M, oshape: OShape,
       side_face: SideFace) -> R;
  
  fn intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(&self,
       mon: M,
       q: &VectorMonomial<M>,
       oshape: OShape,
       side_face: SideFace) -> R;
 

  fn intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side<P:Polynomial<M>>(&self,
       p: &P,
       q: &VectorMonomial<M>,
       oshape: OShape,
       side_face: SideFace) -> R;

}


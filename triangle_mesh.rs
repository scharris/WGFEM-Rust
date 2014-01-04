use common::*;
use monomial::{Monomial};
use polynomial::Polynomial;
use vector_monomial::VectorMonomial;
use mesh::*;
use quadrature::*;
use storage_by_ints::StorageByInts2;

use std::vec;
use std::num::{hypot, min, max, abs};
use std::iter::{Iterator};
use std::cmp::{Less, Greater};


// Reference triangle type representing the oriented shapes of the mesh.
pub struct RefTri {

  v01: Vec, // displacement vector from vertex 0 to vertex 1
  v02: Vec, // displacement vector from vertex 0 to vertex 2
 
  // There can be multiple side faces between vertexes, the numbers of which are recorded here for each vertex pair.
  nums_side_faces_between_vertexes: (u8,u8,u8), // v0 <-> v1, v1 <-> v2, v2 <-> v0
  
  num_side_faces: uint,
  outward_normals_by_side_face: ~[Vec],
  dep_dims_by_side_face: ~[Dim],

  diameter_inv: R,

}

// A single finite element in the final mesh.
pub struct ElTri {

  oshape: OShape, // key identifying a RefTri

  v0: Point,

}

// The mesh.
pub struct TriMesh<Mon> {

  // finite elements, indexed by finite element number
  fes: ~[ElTri],

  // reference triangles representing the oriented shapes of the mesh, indexed by oriented shape number
  oshapes: ~[RefTri],

  // non-boundary side numbers by fe and side face
  nbsidenums_by_fe_face: StorageByInts2<Option<NBSideNum>>,

  // inclusions of non-boundary sides in fe's, indexed by non-boundary side number
  nbsideincls_by_nbsidenum: ~[NBSideInclusions],

  // number of finite elements
  num_fes: uint,

  // number of non-boundary sides
  num_nb_sides: uint,

  // number of boundary sides
  num_b_sides: uint,

  // maximum number of side faces for any oriented shape
  max_num_shape_sides: uint,

  // number of oriented shapes
  num_oshapes: uint,

  // integration support members
  integration_rel_err: R,
  integration_abs_err: R,

  // element -> Tag maps
  phys_reg_tags_by_fenum: Option<~[Tag]>,
  geom_ent_tags_by_fenum: Option<~[Tag]>,

}

pub type Point = (R,R);
pub type Vec = (R,R);

// Tag for a mesh geometric region or physical entity which may be attributed to one or more elements.
#[deriving(Eq, TotalEq, Ord, TotalOrd, Clone)]
pub struct Tag(int);


impl<Mon:Monomial> Mesh<Mon> for TriMesh<Mon> {

  #[inline]
  fn num_fes(&self) -> uint {
    self.num_fes
  }

  #[inline]
  fn num_nb_sides(&self) -> uint {
    self.num_nb_sides
  }

  #[inline]
  fn num_oriented_element_shapes(&self) -> uint {
    self.num_oshapes
  }

  #[inline]
  fn oriented_shape_for_fe(&self, fe: FENum) -> OShape {
    self.fes[*fe].oshape
  }

  #[inline]
  fn num_side_faces_for_oshape(&self, os: OShape) -> uint {
    self.oshapes[*os].num_side_faces
  }

  #[inline]
  fn dependent_dim_for_oshape_side(&self, os: OShape, sf: SideFace) -> Dim {
    self.oshapes[*os].dep_dims_by_side_face[*sf]
  }

  #[inline]
  fn fe_inclusions_of_nb_side(&self, nbsn: NBSideNum) -> NBSideInclusions {
    self.nbsideincls_by_nbsidenum[*nbsn]
  } 

  #[inline]
  // Return non-boundary side number of the indicated fe relative side.
  fn nb_side_num_for_fe_side(&self, fe: FENum, sf: SideFace) -> NBSideNum {
    self.nbsidenums_by_fe_face.get(*fe, *sf).unwrap()
  }

  #[inline]
  fn is_boundary_side(&self, fe: FENum, sf: SideFace) -> bool {
    self.nbsidenums_by_fe_face.get(*fe, *sf).is_none()
  }

  #[inline]
  fn num_boundary_sides(&self) -> uint {
    self.num_b_sides
  }
  
  fn boundary_fes_by_oshape_side(&self) -> ~[~[~[FENum]]] { // oshape, side face -> fes
    let mut b_fes = vec::from_fn(self.num_oshapes, 
                                 |os| vec::from_elem(self.num_side_faces_for_oshape(OShape(os)), ~[]));

    for fe in range(0, self.num_fes) {
      let fe = FENum(fe);
      let os = self.oriented_shape_for_fe(fe);
      for sf in range(0, self.num_side_faces_for_oshape(os)) {
        if self.is_boundary_side(fe, SideFace(sf)) {
          b_fes[*os][sf].push(fe);
        }
      }
    }
    
    b_fes
  }
  
  #[inline]
  fn shape_diameter_inv(&self, os: OShape) -> R {
    self.oshapes[*os].diameter_inv
  }

  #[inline]
  fn max_fe_diameter(&self) -> R {
    self.oshapes.iter().map(|ref_tri|  1./ref_tri.diameter_inv).max().unwrap()
  }

  #[inline]
  fn num_nb_sides_for_fe(&self, fe: FENum) -> uint {
    range(0, self.num_side_faces_for_oshape(self.oriented_shape_for_fe(fe)))
      .count(|sf| !self.is_boundary_side(fe, SideFace(sf))) 
  }

  #[inline]
  fn max_num_shape_sides(&self) -> uint {
    self.max_num_shape_sides
  }

  // integration functions
  
  fn intg_global_fn_on_fe_int
     ( &self,
       f:  |&[R]| -> R,
       fe: FENum )
     -> R
  {
    // Sort the points so we can divide the triangle into regions between sloped lines diverging from a point above
    // and below and a vertical line on one side.
    let (p0, p1, p2) = self.sorted_vertexes(self.ref_tri_for_fe(fe), self.fes[*fe].v0);
    
    let (slope_01, slope_02) = (slope_between(p0, p1), slope_between(p0, p2));

    if p0.n0() == p1.n0() { // Points 0 and 1 are on a vertical line on the left, point 2 on the right.
     /*  p1
      *  |\
      *  | \ p2
      *  | /
      *  |/
      *  p0
      *
      * Integrate over the area bounded above and below by the lines between points 0 and 2 and 1 and 2, and
      * horizontally between the vertical left side formed by points 0 and 1, and point 2 on the right where
      * the upper and lower bounding lines meet.
      */
      self.intg_fn_btw_pt_and_vert_seg(f, p2, slope_02, slope_between(p1, p2), p0.n0())
    } else {
     /* Points 0 and 1 do not form a vertical line.
      *      p1            p2
      *     /|\           /|
      *    / | \         / |
      * p0/__|__\p2   p0/__|p1  (bottom lines not necessarily horizontal)
      *
      * Integrate between points 0 and 1, and between points 1 and 2 if points 1 and 2 don't lie on a vertical line.
      */
      let left_seg = self.intg_fn_btw_pt_and_vert_seg(|x| f(x), p0, slope_01, slope_02, p1.n0());
      let right_seg = {
        if p1.n0() == p2.n0() { 0 as R } // vertical right side, no second segment
        else {
          self.intg_fn_btw_pt_and_vert_seg(f, p2, slope_02, slope_between(p1, p2), p1.n0())
        }
      };
      left_seg + right_seg
    }
  }

  // Integrate a global function against an interior-relative monomial over a finite element interior.
  #[inline]
  fn intg_global_fn_x_facerel_mon_on_fe_int
     ( &self,
       f:   |&[R]| -> R,
       mon: Mon,
       fe:  FENum )
     -> R
  {
    let int_origin_arr = { let int_origin = self.fes[*fe].v0; [int_origin.n0(), int_origin.n1()] };
    let integrand = |x: &[R]| { f(x) * mon.value_at_for_origin(x, int_origin_arr) };
    self.intg_global_fn_on_fe_int(integrand, fe) 
  }

  fn intg_global_fn_x_facerel_mon_on_fe_side
     ( &self,
       f:   |&[R]| -> R,
       mon: Mon,
       fe:  FENum,
       sf:  SideFace )
     -> R
  {
    let ref_tri = self.ref_tri_for_fe(fe);
    let v0 = self.fes[*fe].v0;
    let (v1, v2) = (vsum(v0, ref_tri.v01), vsum(v0, ref_tri.v02));
    
    // We want to compute int_0^1 f(s(t)) mon(s(t)-o) |s'(t)| dt, where s:[0,1] -> R^2 is a bijection traversing the
    // side face smoothly with s(0) and s(1) being the side endpoints, and o is the local origin for the side, which
    // is the side's midpoint.
    let (a,b) = side_face_endpoint_pair(sf, v0,v1,v2, ref_tri.nums_side_faces_between_vertexes, VertexOrder);
    let o = { let mp = midpt(a, b); [mp.n0(), mp.n1()] }; // The local origin of each side face is the midpoint.
    let side_len = dist(a,b);
      
    let mut s_t = [0.,0.];
    let integrand = |t: &[R]| { // compute f(s(t)) mon(s(t)-o) |s'(t)|
      let t = t[0];
      s_t[0] = a.n0() + (b.n0() - a.n0())*t;
      s_t[1] = a.n1() + (b.n1() - a.n1())*t;
      f(s_t) * mon.value_at_for_origin(s_t, o) * side_len  // side_len = |s'(t)|
    };

    space_adaptive_quadrature(&integrand,
                              [0. as R], [1. as R],
                              self.integration_rel_err, self.integration_abs_err)
  }
  
  fn intg_mixed_global_and_facerel_fn_on_fe_int
     ( &self,
       f:  |&[R], &[R]| -> R,
       fe: FENum )
     -> R
  {
    let int_origin = self.fes[*fe].v0;

    let mut x_rel = [0. as R, 0. as R]; 
    let integrand = |x: &[R]| {
      x_rel[0] = x[0] - int_origin.n0();
      x_rel[1] = x[1] - int_origin.n1();
      f(x, x_rel)
    };

    self.intg_global_fn_on_fe_int(integrand, fe)
  }

  #[inline]
  fn intg_facerel_poly_on_oshape_int <P: Polynomial<Mon>>
     ( &self,
       p:  &P,
       os: OShape)
     -> R
  {
    self.intg_facerel_poly_fn_on_oshape_int(|x| p.value_at(x), p.max_var_deg(), os)
  }
  
  #[inline]
  fn intg_facerel_mon_on_oshape_int
     ( &self,
       mon: Mon,
       os:  OShape)
     -> R
  {
    self.intg_facerel_poly_fn_on_oshape_int(|x| mon.value_at(x), mon.max_var_deg(), os)
  }


  #[inline]
  fn intg_facerel_poly_x_facerel_poly_on_oshape_side <P: Polynomial<Mon>>
     ( &self,
       p1: &P,
       p2: &P,
       os: OShape,
       sf: SideFace )
     -> R
  {
    self.intg_facerel_poly_fn_on_oshape_side(|x| p1.value_at(x) * p2.value_at(x), Deg(*p1.deg() + *p2.deg()), os, sf)
  }
  
  #[inline]
  fn intg_facerel_mon_on_oshape_side
     ( &self,
       mon: Mon,
       os:  OShape,
       sf:  SideFace )
     -> R
  {
    self.intg_facerel_poly_fn_on_oshape_side(|x| mon.value_at(x), mon.deg(), os, sf)
  }


  #[inline]
  fn intg_facerel_mon_x_facerel_poly_on_oshape_side <P: Polynomial<Mon>>
     ( &self,
       mon: Mon,
       p:   &P,
       os:  OShape,
       sf:  SideFace)
     -> R
  {
    self.intg_facerel_poly_fn_on_oshape_side(|x| mon.value_at(x) * p.value_at(x), Deg(*mon.deg() + *p.deg()), os, sf)
  }

  fn intg_intrel_mon_x_siderel_mon_on_oshape_side
     ( &self,
       int_mon:  Mon,
       side_mon: Mon,
       os:       OShape,
       sf:       SideFace )
     -> R
  {
    let ref_tri = &self.oshapes[*os];
    let (v0, v1, v2) = ((0.,0.), ref_tri.v01, ref_tri.v02);
    let (a,b) = side_face_endpoint_pair(sf, v0,v1,v2, ref_tri.nums_side_faces_between_vertexes, VertexOrder);
    let sf_o = midpt(a, b); // local origin of the side face relative to v0

    let p = |x_sf: &[R]| { // side relative point (side midpoint origin)
      let x_int = [x_sf[0] + sf_o.n0(), x_sf[1] + sf_o.n1()];
      int_mon.value_at(x_int) + side_mon.value_at(x_sf)
    };
    
    let p_deg = Deg(*int_mon.deg() + *side_mon.deg());

    self.intg_facerel_poly_fn_on_oshape_side(p, p_deg, os, sf)
  }
  
  fn intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side
     ( &self,
       mon: Mon,
       q:   &VectorMonomial<Mon>,
       os:  OShape,
       sf:  SideFace)
     -> R
  {
    let ref_tri = &self.oshapes[*os];
    let (v0, v1, v2) = ((0.,0.), ref_tri.v01, ref_tri.v02);
    let (a,b) = side_face_endpoint_pair(sf, v0,v1,v2, ref_tri.nums_side_faces_between_vertexes, VertexOrder);
    let sf_o = midpt(a, b); // local origin of the side face relative to v0

    let onormal = { let n = ref_tri.outward_normals_by_side_face[*sf]; [n.n0(), n.n1()] };

    let p = |x_sf: &[R]| { // side relative point (side midpoint origin)
      let x_int = [x_sf[0] + sf_o.n0(), x_sf[1] + sf_o.n1()];
      mon.value_at(x_sf) + q.dot_at(onormal, x_int)
    };

    let ub_p_deg = Deg(*mon.deg() + *q.mon.deg());

    self.intg_facerel_poly_fn_on_oshape_side(p, ub_p_deg, os, sf)
  }
 
} // Mesh<Mon> impl


impl<Mon:Monomial> TriMesh<Mon> {

  pub fn new
         ( fes: ~[ElTri],
           oshapes: ~[RefTri],
           nbsidenums_by_fe_face: StorageByInts2<Option<NBSideNum>>,
           nbsideincls_by_nbsidenum: ~[NBSideInclusions],
           num_b_sides: uint,
           integration_rel_err: R,
           integration_abs_err: R,
           phys_reg_tags_by_fenum: Option<~[Tag]>,
           geom_ent_tags_by_fenum: Option<~[Tag]> )
         -> TriMesh<Mon>
  {
    let num_fes = fes.len(); 
    let num_oshapes = oshapes.len();
    let num_nb_sides = nbsideincls_by_nbsidenum.len();
    let max_num_shape_sides = oshapes.iter().map(|ref_tri| ref_tri.num_side_faces).max().unwrap();

    TriMesh {
      fes: fes,
      oshapes: oshapes,
      nbsidenums_by_fe_face: nbsidenums_by_fe_face,
      nbsideincls_by_nbsidenum: nbsideincls_by_nbsidenum,
      num_fes: num_fes,
      num_nb_sides: num_nb_sides,
      num_b_sides: num_b_sides,
      max_num_shape_sides: max_num_shape_sides,
      num_oshapes: num_oshapes,
      integration_rel_err: integration_rel_err,
      integration_abs_err: integration_abs_err,
      phys_reg_tags_by_fenum: phys_reg_tags_by_fenum,
      geom_ent_tags_by_fenum: geom_ent_tags_by_fenum,
    }
  }
    

  // Reference triangle for a finite element.
  fn ref_tri_for_fe<'a>(&'a self, fe: FENum) -> &'a RefTri
  {
    &self.oshapes[*self.fes[*fe].oshape]
  }


 /* Integrate a function over the triangular region bounded on two sides by two non-vertical lines of
  * indicated slopes which meet at a point q, and by the indicated vertical line as the remaining side.
  */
  fn intg_fn_btw_pt_and_vert_seg
     ( &self,
       f:           |&[R]| -> R,
       q:           Point,
       slope_1:     R,
       slope_2:     R,
       vert_line_x: R )
     -> R
  {
    let (xminT, xmaxT) = (min(q.n0(), vert_line_x), max(q.n0(), vert_line_x));
    let w = xmaxT - xminT; // width of triangular integration region

    // Method
    // We will pull back the integration over the original triangular section T to an integration
    // over the unit square, by change of variables via the bijection
    // 
    //   t:[0,1]^2 -> T:  t(x,y) = (xminT + x w,  q_2 + m1 (xminT + x w - q_1) + y(m2 - m1)(xminT + x w - q_1))
    // 
    // Here m1 and m2 are the slopes of the non-vertical bounding lines (m1 != m2).
    // The determinant of the derivative matrix Dt(x,y) is
    //                     |          w                           0                |
    //  det Dt(x,y)| = det |                                                       |
    //                     | m1 w + y(m2 - m1) w      (m2 - m1)(xminT + x w - q_1) |
    //               = w (m2 - m1) (xminT + x w - q_1).
    //  Now by applying change of variables in the integral of f over the T via the mapping t,
    //  we have
    //    int_T f = int_0^1 int_0^1 f(t(x,y)) |det Dt(x,y)| dy dx
    //            = int_0^1 int_0^1 f(t(x,y)) |w (m2 - m1) (xminT + x w - q_1)| dy dx
    
    let slopediff = slope_2 - slope_1;
    let w_slopediff = w * slopediff;
    
    let mut t = [0.,0.];
    let integrand = |s: &[R]| { // unit square point
      let xT = xminT + s[0] * w;     // triangle x
      let xT_minus_qx = xT - q.n0(); // relative triangle x
        
      // Construct the triangle point t(x,y).
      t[0] = xT;
      t[1] = q.n1() + slope_1 * xT_minus_qx + s[1] * slopediff * xT_minus_qx;
      let det_Dt = w_slopediff * xT_minus_qx;
      f(t) * abs(det_Dt)
    };
    
    space_adaptive_quadrature(&integrand,
                              [0 as R, 0 as R], [1 as R, 1 as R],
                              self.integration_rel_err, self.integration_abs_err)
  }

  fn intg_facerel_poly_fn_on_oshape_int
     ( &self,
       p:                   |&[R]| -> R,
       ub_max_var_deg_in_p: Deg,
       os:                  OShape )
     -> R
  {
    // Sort the points so we can divide the triangle into regions between sloped lines diverging from a point above
    // and below and a vertical line on one side.
    let (p0, p1, p2) = self.sorted_vertexes(&self.oshapes[*os], (0.,0.)); // vertexes relative to interior origin, which is vertex 0.
    let (slope_01, slope_02) = (slope_between(p0, p1), slope_between(p0, p2));
   
    if p0.n0() == p1.n0() { // Points 0 and 1 are on a vertical line on the left, point 2 on the right.
     /*  p1
      *  |\
      *  | \ p2
      *  | /
      *  |/
      *  p0
      */
      let slope_12 = slope_between(p1, p2);
      self.intg_poly_fn_btw_pt_and_vert_seg(p, ub_max_var_deg_in_p, p2, slope_02, slope_12, p0.n0())
    } else {
     /*  Points 0 and 1 do not form a vertical line.      
      *      p1            p2
      *     /|\           /|
      *    / | \         / |
      * p0/__|__\p2   p0/__|p1  (bottom lines not necessarily horizontal)
      */
      let left_seg = self.intg_poly_fn_btw_pt_and_vert_seg(|x| p(x), ub_max_var_deg_in_p, p0, slope_01, slope_02, p1.n0());
      let right_seg = {
        if p1.n0() == p2.n0() { 0 as R } // vertical right side, no second segment
        else
        {
          let slope_12 = slope_between(p1, p2);
          self.intg_poly_fn_btw_pt_and_vert_seg(p, ub_max_var_deg_in_p, p2, slope_02, slope_12, p1.n0())
        }
      };
      left_seg + right_seg
    }
  }


 /* Integrate a polynomial as a global function (without its own coordates origin) over the triangular region bounded on
  * two sides by two non-vertical lines of indicated slopes which meet at a point q, and by the indicated vertical line
  * as the remaining side.
  */
  fn intg_poly_fn_btw_pt_and_vert_seg
     ( &self, 
       p:                   |&[R]| -> R,
       ub_max_var_deg_in_p: Deg,
       q:                   Point,
       slope_1:             R,
       slope_2:             R,
       vert_line_x:         R )
     -> R
  {
    let (xminT, xmaxT) = (min(q.n0(), vert_line_x), max(q.n0(), vert_line_x));
    let w = xmaxT - xminT; // width of triangular integration region

    let slopediff = slope_2 - slope_1;
    let w_slopediff = w * slopediff;
    
    // The integrand after change of variables to pull back the integration to the unit square is:
    //   (x,y) -> p(t(x,y)) |w (m2 - m1) (xminT + x w - q_1)|
    //   where t:[0,1]^2 -> T:  t(x,y) = (xminT + x w,  q_2 + m1 (xminT + x w - q_1) + y(m2 - m1)(xminT + x w - q_1)).
    // See comments in above function about the change of variables being used.
    let mut t = [0.,0.];
    let unit_sq_integrand = |xS: R, yS: R| { // unit square point
      let xT = xminT + xS * w;       // triangle x
      let xT_minus_qx = xT - q.n0(); // relative triangle x
      // Construct the triangle point t(x,y).
      t[0] = xT;
      t[1] = q.n1() + slope_1 * xT_minus_qx + yS * slopediff * xT_minus_qx;
      let det_Dt = w_slopediff * xT_minus_qx;
      p(t) * abs(det_Dt)
    };
    
    // The integrand is a polynomial because the RHS factor under the absolute value will not change sign
    // for x in [0,1]: its non-constant factor will range either from 0 to w if q_1 = xminT, or from -w to
    // 0 if q_1 = xmaxT.
    // 
    // For the integral to be exact using repeated Gaussian quadrature, we need
    //   k <= 2n - 1,
    // where 
    //   n is the number of weight values for each axis (2n weights, n^2 total points), and
    //   k is the maximum individual variable degree in the integrand polynomial.
    //   [See e.g. http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF]
    // The maximum variable degree in (x,y) -> p(t(x,y)) is no more than the maximum variable degree of p itself,
    // because t's maximum variable degree is 1 in both components. Taking the right hand absolute value factor
    // into account, which has a maximum variable degree of 1, then an upper bound k on the maximum variable
    // degree for the complete integrand is one more than the maximum variable degree in p.
    let k = (1 + *ub_max_var_deg_in_p) as uint;
    
    // Choose the smallest number of weights per axis n such that k <= 2n - 1.
    let n = if k % 2 != 0 { (k+1)/2 } else { (k+1)/2 + 1 };
    
    gaussian_quadrature_2D_rect(n, &unit_sq_integrand, 0 as R, 0 as R, 1 as R, 1 as R)
  }

  fn intg_facerel_poly_fn_on_oshape_side
     ( &self,
       p:        |&[R]| -> R,
       ub_p_deg: Deg,
       os:       OShape,
       sf:       SideFace )
     -> R
  {
    // We want to compute 
    //   int_0^1 p(s(t)-o) |s'(t)| dt,
    // where s:[0,1] -> R^2 is a bijection traversing the side face smoothly, s(0) and s(1) are the side endpoints,
    // and o is the local origin for the side which is the side's midpoint.

    let ref_tri = &self.oshapes[*os];
    let (v0, v1, v2) = ((0.,0.), ref_tri.v01, ref_tri.v02);

    let (a,b) = side_face_endpoint_pair(sf, v0,v1,v2, ref_tri.nums_side_faces_between_vertexes, VertexOrder);
    let o = midpt(a, b); // The local origin of each side face is the midpoint.
    let side_len = dist(a,b);
     
    let integrand = |t: R| { // compute p(s(t)-o) |s'(t)|
      let s_t_minus_o = [a.n0() + (b.n0() - a.n0())*t - o.n0(),
                         a.n1() + (b.n1() - a.n1())*t - o.n1()];
      p(s_t_minus_o) * side_len  // side_len = |s'(t)|
    };

    // The integrand degree will be at most that of p itself, because the path s has a polynomial of degree at most 1
    // in each of its output components.
    let k = *ub_p_deg as uint;
    
    // For an exact result, we need a number of quadrature points n such that k <= 2n - 1.
    let n = if k % 2 != 0 { (k+1)/2 } else { (k+1)/2 + 1 };
    
    gaussian_quadrature(n, &integrand, 0 as R, 1 as R)
  }

  #[inline]
  fn sorted_vertexes
     ( &self,
       ref_tri: &RefTri,
       v0:      Point )
     -> (Point,Point,Point)
  {
    let sorted_pts = {
      let mut pts = [v0, vsum(v0, ref_tri.v01), vsum(v0, ref_tri.v02)];
      pts.sort_by(|p,q| if p.n0() < q.n0() || p.n0() == q.n0() && p.n1() < q.n1() { Less } else { Greater }); // (we know points are not eq)
      pts
    };
    (sorted_pts[0], sorted_pts[1], sorted_pts[2])
  }

} // TriMesh<Mon> impl


impl RefTri {

  pub fn new
         ( v0: Point,
           v1: Point, 
           v2: Point,
           scale: R, 
           nums_side_faces_btw_verts: (u8,u8,u8) )
         -> RefTri
  {
    let ((v0_x,v0_y), (v1_x,v1_y), (v2_x,v2_y)) = (v0,v1,v2);
    let scaled_v01 = (scale*(v1_x - v0_x), scale*(v1_y - v0_y));
    let scaled_v02 = (scale*(v2_x - v0_x), scale*(v2_y - v0_y));
    let norm_scaled_v12 = scale * dist(v1, v2);
    let normals = RefTri::side_face_outward_normals(v0, v1, v2, nums_side_faces_btw_verts);
    let dep_dims_by_side_face = RefTri::side_face_dep_dims(v0,v1,v2, nums_side_faces_btw_verts);
    let diameter_inv = 1./max(norm(scaled_v01), max(norm(scaled_v02), norm_scaled_v12));
    let num_side_faces = {
      let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
      (sfs_v01 + sfs_v12 + sfs_v20) as uint
    };
    RefTri{ v01: scaled_v01,
            v02: scaled_v02,
            nums_side_faces_between_vertexes: nums_side_faces_btw_verts,
            outward_normals_by_side_face: normals,
            dep_dims_by_side_face: dep_dims_by_side_face,
            diameter_inv: diameter_inv,
            num_side_faces: num_side_faces }
  }

 /*
  * Outward Normals for Side Faces
  * The outward normal for an inter-vertex vector w (extended to R^3 via 0 third component) is
  *   n = (w x (0,0,cc)) / |w|
  *     = cc (w_2, -w_1, 0) / |w|. // (using one based indexing on w)
  * where cc = 1 if w is part of a clockwise traversal of the triangle, and -1 otherwise.
  * For any pair of successive inter-vertex vectors v and w, cc can be computed as:
  *   cc = sgn(z(v x w)), where z is projection of component 3.
  */
  fn side_face_outward_normals
     ( v0: Point,
       v1: Point,
       v2: Point,
       nums_side_faces_btw_verts: (u8,u8,u8) )
     -> ~[Vec]
  {
    let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
    // inter-vertex vectors
    let (v01, v12, v20) = (vdiff(v1,v0), vdiff(v2,v1), vdiff(v0,v2));
    let ivvs = [(v01, sfs_v01),
                (v12, sfs_v12),
                (v20, sfs_v20)];
    let cc = if oriented_counterclockwise(v01, v12) { 1. as R } else { -1. as R };
    
    let mut normals = vec::with_capacity((sfs_v01 + sfs_v12 + sfs_v20) as uint);
    for &((ivv_x, ivv_y), sfs) in ivvs.iter() {
      let len = hypot(ivv_x, ivv_y);
      let n = (cc * ivv_y/len, -cc * ivv_x/len);
      (sfs as uint).times(|| normals.push(n));
    }
    normals
  }
  
  fn side_face_dep_dims
     ( v0: Point,
       v1: Point,
       v2: Point,
       nums_side_faces_btw_verts: (u8,u8,u8) )
     -> ~[Dim]
  {
    let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
    let vert_pairs = [((v0,v1), sfs_v01),
                      ((v1,v2), sfs_v12),
                      ((v2,v0), sfs_v20)];
    
    let mut dep_dims = vec::with_capacity((sfs_v01 + sfs_v12 + sfs_v20) as uint);

    for &((va,vb), sfs) in vert_pairs.iter() {
      let dep_dim = RefTri::side_face_dep_dim(va, vb);
      (sfs as uint).times(|| dep_dims.push(dep_dim));
    }
    dep_dims
  }

  // Return a dependent dimension for the given side endpoints.
  fn side_face_dep_dim(va: Point, vb: Point) -> Dim {
    if abs(vb.n0() - va.n0()) >= abs(vb.n1() - va.n1()) { Dim(1) } else { Dim(0) }
  }

} // RefTri impl


// standalone auxiliary functions and related types

pub enum SideEndpointsOrdering {
  LesserEndpointsFirst,
  VertexOrder,
}

fn side_face_endpoint_pair
   ( sf: SideFace,
     v0: Point,
     v1: Point,
     v2: Point,
     nums_side_faces_btw_verts: (u8,u8,u8),
     endpoints_ordering: SideEndpointsOrdering )
   -> (Point, Point)
{
  let mut cur_sf = 0u; 
  
  let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
  let side_faces = [(v0, v1, sfs_v01),
                    (v1, v2, sfs_v12),
                    (v2, v0, sfs_v20)];
  
  for &(va, vb, num_sfs_btw_va_vb) in side_faces.iter() {
    match num_sfs_btw_va_vb {
      1 => {
        if cur_sf == *sf { return mk_endpoint_pair(va, vb, endpoints_ordering) }
        cur_sf += 1;
      }
      2 => {
        let midpt = midpt(va, vb);
        if cur_sf == *sf { return mk_endpoint_pair(va, midpt, endpoints_ordering); }
        cur_sf += 1;
        if cur_sf == *sf { return mk_endpoint_pair(midpt, vb, endpoints_ordering); }
        cur_sf += 1;
      }
      _ =>  fail!("Only 1 or 2 faces between triangle vertexes are currently supported.")
    }
  }
  
  fail!("Side face $sf not found.");
}

#[inline]
pub fn mk_endpoint_pair
       ( pt1:      Point,
         pt2:      Point,
         ordering: SideEndpointsOrdering )
       -> (Point, Point)
{
  match ordering {
    LesserEndpointsFirst => if pt1 < pt2 { (pt1, pt2) } else { (pt2, pt1) },
    VertexOrder => (pt1, pt2)
  }
}

#[inline]
fn norm((v_x, v_y): Vec) -> R { hypot(v_x, v_y) }

#[inline]
fn dist(p: Point, q: Point) -> R { hypot(q.n0() - p.n0(), q.n1() - p.n1()) }

#[inline]
fn slope_between(p: Point, q: Point) -> R { (p.n1() - q.n1())/(p.n0() - q.n0()) }

#[inline]
pub fn midpt(p: Point, q: Point) -> Point { (0.5*(p.n0() + q.n0()), 0.5*(p.n1() + q.n1())) }

#[inline]
fn vsum(p: Vec, q: Vec) -> Vec { (p.n0() + q.n0(), p.n1() + q.n1()) }

#[inline]
fn vdiff(p: Vec, q: Vec) -> Vec { (p.n0() - q.n0(), p.n1() - q.n1()) }
  
// Determine whether a counterclockwise rotation not exceeding 1/2 turn can transform the first vector to a positive
// multiple of the second.
#[inline]
fn oriented_counterclockwise((u_x,u_y): Vec, (v_x,v_y): Vec) -> bool { u_x*v_y - u_y*v_x > 0. as R }


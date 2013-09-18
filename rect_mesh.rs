extern mod extra;
use extra::treemap::TreeMap;
use std::ptr;
use std::vec;
use std::num::{sqrt, abs};
use std::iter::range_inclusive;
use common::*;
use monomial::{Monomial, Mon1d, Mon2d, Mon3d, Mon4d};
use polynomial::Polynomial;
use vector_monomial::VectorMonomial;
use mesh::*;
use quadrature::*;

mod common;
mod monomial;
mod polynomial;
mod vector_monomial;
mod mesh;
mod quadrature;

static DEFAULT_INTEGRATION_REL_ERR: R = 1e-12;
static DEFAULT_INTEGRATION_ABS_ERR: R = 1e-12;

// MeshCoord represents a single integer-valued mesh coordinates component,
// e.g. column or row or stack, etc.
#[deriving(Eq, TotalEq, Ord, TotalOrd, Clone)]
pub struct MeshCoord(uint);

pub struct NBSideGeom {
  perp_axis: Dim,
  mesh_coords: ~[MeshCoord]
}

pub struct RectMesh<M> {

  // The dimension of the Euclidiean space containing the mesh. 
  space_dim: Dim,

  // Mesh coordinate ranges in R^d defining the boundaries of the mesh.
  min_bounds: ~[R],
  max_bounds: ~[R],

  // Logical dimensions of the mesh in integer mesh axis coordinates,
  // with directions corresponding to the coordinate axes (cols, rows,...).
  mesh_ldims: ~[MeshCoord],

  // Actual dimensions of any single finite element, the displacement vector from the
  // minimum coordinates corner to the maximum coordinates corner.
  fe_dims: ~[R],

  // Cumulative products of the logical dimensions of the main finite element mesh. The r^th 
  // component is the product of the logical dimensions up to and including dimension r.
  cumprods_mesh_ldims: ~[uint],

  // Cumulative products of the logical dimensions of sides meshes, by side perpendicular axis.
  // Component a holds the vector of cumulative products of logical dimensions of the mesh of sides
  // perpendicular to axis a. A side mesh's cumulative products at component r is the product of
  // logical dimensions of the mesh up to and including dimension r.
  cumprods_nb_side_mesh_ldims_by_perp_axis: ~[~[uint]],

  // An assignment of unique side numbers ranges to groups of non-boundary sides, where the sides
  // are ordered and grouped by their perpendicular axis number in increasing order.
  first_nb_side_nums_by_perp_axis: ~[NBSideNum],

  // The number of finite elements in the mesh.
  num_fes: uint,
  // The number of non-boundary sides in the mesh.
  num_nb_sides: uint,
  // The number of side faces (boundary or non-boundary) for each finite element.
  num_side_faces_per_fe: uint,

  // The diameter (diagonal length) of the element rectangles in the mesh.
  rect_diameter: R,
  // 1/rect_diameter  
  rect_diameter_inv: R,

  // The tolerated relative and absolute erros for numerical integration.
  integration_rel_err: R,
  integration_abs_err: R
}


fn new_impl<M:Monomial>(min_bounds: ~[R],
                        max_bounds: ~[R],
                        mesh_ldims: ~[MeshCoord],
                        integration_rel_err: R,
                        integration_abs_err: R) -> ~RectMesh<M> {

  let space_dim = Monomial::domain_dim(None::<M>);
  assert!(min_bounds.len() == *space_dim);
  assert!(max_bounds.len() == *space_dim);
  assert!(mesh_ldims.len() == *space_dim);
  
  let fe_dims: ~[R] =
    vec::from_fn(*space_dim, |r| {
      let bounds_diff = max_bounds[r] - min_bounds[r];
      let ldim_r = *mesh_ldims[r];
      assert!(bounds_diff > 0 as R);
      assert!(ldim_r > 0);
      bounds_diff/(ldim_r as R)
    });

  let fe_dims_with_drops: ~[~[R]] =
    vec::from_fn(*space_dim, |r| {
      if r != *space_dim-1 { fe_dims.slice_to(r) + fe_dims.slice_from(r+1) }
      else { fe_dims.slice_to(r).to_owned() }
    });
 
  let cumprods_mesh_ldims: ~[uint] =
    mesh_ldims.iter().scan(1, |prod, &ldim| {
      *prod *= *ldim;
      Some(*prod)
    }).to_owned_vec();

  let cumprods_nb_side_mesh_ldims_by_perp_axis: ~[~[uint]] =
    vec::from_fn(*space_dim, |perp_axis| {
      vec::from_fn(*space_dim, |prods_top_dim| {
        let mut prod = 1u;
        for r in range_inclusive(0, prods_top_dim) {
          prod *= if r != perp_axis { *mesh_ldims[r] } else { *mesh_ldims[r]-1 }
        }
        prod
      })
    });

  let num_fes = *cumprods_mesh_ldims.last();

  let nb_side_counts_by_perp_axis = cumprods_nb_side_mesh_ldims_by_perp_axis.iter()
                                      .map(|cumprods| *cumprods.last())
                                      .to_owned_vec();

  let num_nb_sides = nb_side_counts_by_perp_axis.iter().fold(0u, |sum, &x| sum + x);

  let first_nb_side_nums_by_perp_axis: ~[NBSideNum] = {
    ~[NBSideNum(0u)] +
    nb_side_counts_by_perp_axis.init()
      .iter().scan(0, |sum, &axis_nb_sides| { *sum += axis_nb_sides; Some(NBSideNum(*sum)) })
      .to_owned_vec()
  };

  let rect_diameter = sqrt(fe_dims.iter().fold(0 as R, |sum_sq_dims, &fe_dim| sum_sq_dims + fe_dim*fe_dim));

  ~RectMesh {
    space_dim: space_dim,
    min_bounds: min_bounds,
    max_bounds: max_bounds,
    mesh_ldims: mesh_ldims,
    fe_dims: fe_dims,
    cumprods_mesh_ldims: cumprods_mesh_ldims,
    cumprods_nb_side_mesh_ldims_by_perp_axis: cumprods_nb_side_mesh_ldims_by_perp_axis,
    first_nb_side_nums_by_perp_axis: first_nb_side_nums_by_perp_axis,
    num_fes: num_fes,
    num_nb_sides: num_nb_sides,
    num_side_faces_per_fe: 2 * *space_dim,
    rect_diameter: rect_diameter,
    rect_diameter_inv: 1./rect_diameter,
    integration_rel_err: integration_rel_err,
    integration_abs_err: integration_abs_err
  }
}


impl<M:Monomial> RectMesh<M> {
  
  pub fn new(min_bounds: ~[R],
             max_bounds: ~[R],
             mesh_ldims: ~[MeshCoord]) -> ~RectMesh<M> {
    new_impl(min_bounds, max_bounds, mesh_ldims,
             DEFAULT_INTEGRATION_REL_ERR, DEFAULT_INTEGRATION_ABS_ERR)
  }

  pub fn new_with_err_tols(min_bounds: ~[R],
                           max_bounds: ~[R],
                           mesh_ldims: ~[MeshCoord],
                           integration_rel_err: R,
                           integration_abs_err: R) -> ~RectMesh<M> {
      new_impl(min_bounds, max_bounds, mesh_ldims,
               integration_rel_err, integration_abs_err)
  }


  // side-related functions

  // Finds the perpendicular axis for a given non-boundary side in the mesh.
  #[inline]
  fn perp_axis_for_nb_side(&self, n: NBSideNum) -> Dim {
    assert!(*n < self.num_nb_sides);
    let mut r = *self.space_dim-1;
    loop {
      if self.first_nb_side_nums_by_perp_axis[r] <= n { return Dim(r) }
      else if r == 0 { fail!("cannot find perpendicular axis for non-boundary side"); }
      else { r -= 1 };
    }
  }
  
  // Returns the geometric information for a non-boundary side in the mesh, which identifies the
  // perpendicular axis for the side, together with its coordinates in the orientation-specific mesh.
  fn nb_side_geom(&self, n: NBSideNum) -> NBSideGeom {
    /*
     * The r^th logical 0-based coordinate of side n in the mesh of sides having the same orientation is
     *   π(r,n) = ((n − s_a(n)) mod Prod_{i=0..r} k_{a(n),i}) \ (Prod_{i=0..r-1} k_{a(n),i}), (r = 0,...,d-1)
     * where
     *   s_a is the number of the first side in the nb-side enumeration perpendicular to axis a
     *   a(n) is the axis number to which side n is perpendicular
     *   k_{a,i} is component i of the dimensions of the mesh of sides perpendicular to axis a
     * See Rectangular_Meshes.pdf document for the derivation.
     */
    let a = self.perp_axis_for_nb_side(n);
    let orientation_rel_side_num = *n - *self.first_nb_side_nums_by_perp_axis[*a];
    let cumprods_ldims_through = &self.cumprods_nb_side_mesh_ldims_by_perp_axis[*a];
    let sidemesh_coords = vec::from_fn(*self.space_dim, |r| {
      let cumprods_preceeding_ldims = if r == 0 { 1 } else { cumprods_ldims_through[r-1] };
      MeshCoord((orientation_rel_side_num % cumprods_ldims_through[r]) / cumprods_preceeding_ldims)
    });

    NBSideGeom {
      perp_axis: a,
      mesh_coords: sidemesh_coords
    }
  }

  // Converts a side with a given perpendicular axis and orientation-specific side mesh coordinates
  // to a non-boundary side number.
  fn nb_side_with_mesh_coords(&self, coords: &[MeshCoord], perp_axis: Dim) -> NBSideNum {
    /*
     * The enumeration number for a non-boundary side perpendicular to a given axis a, with mesh
     * coordinates (c_1,...,c_d) in its orientation-specific non-boundary side mesh, is
     *   s_{a,#}(c_1,...,c_d) = s_a0 + sum_{r=1..d} { c_r prod_{l=1..r-1} k_{a,l} }
     *                       = s_a0 + c_1 + sum_{r=2..d} { c_r prod_{l=1..r-1} k_{a,l} }
     * where
     *  s_a0 is the enumeration number of the first axis-a perpendicular non-boundary side,
     *  k_{a,l} is the l^th component of the dimensions of the mesh of non-boundary sides
     *    which are perpendicular to axis a.
     */
    let s_a0 = self.first_nb_side_nums_by_perp_axis[*perp_axis];
    let mut sum = *s_a0 + *coords[0];
    for r in range(1, *self.space_dim) {
      sum += *coords[r] * self.cumprods_nb_side_mesh_ldims_by_perp_axis[*perp_axis][r-1];
    }
    NBSideNum(sum)
  }

  // Converts finite element/interior coords in the main mesh to a finite element/interior number.
  #[inline]
  fn fe_with_mesh_coords(&self, coords: &[MeshCoord]) -> FENum {
    /*
     * The finite element (or interior) number for given mesh coordinates (c_1,...,c_d) is
     *   i_#(c_1,...,c_d) = sum_{i=1..d} { c_i prod_{l=1..i-1} k_l }
     *                   = c_1 + sum_{i=2..d} { (c_i - 1) prod_{l=1..i-1} k_l }
     * where k_l is the l^th component of the mesh dimensions.
     * Sum the contributions from the coordinates, each coordinate's value weighted with the
     * cumulative product of the logical mesh dimensions for lesser coordinate dimensions.
     * See Rectangular_Meshes.pdf document for the derivation.
     */
    let mut coord_contrs = *coords[0];
    for r in range(1, *self.space_dim) {
      coord_contrs += *coords[r] * self.cumprods_mesh_ldims[r-1];
    }
    FENum(coord_contrs)
  }

  #[inline]
  fn fe_mesh_coords(&self, fe: FENum) -> ~[MeshCoord] {
    vec::from_fn(*self.space_dim, |r| self.fe_mesh_coord(Dim(r), fe))
  }
  #[inline]
  fn fe_mesh_coord(&self, r: Dim, fe: FENum) -> MeshCoord {
    /*
     * The r^th 0-based mesh coordinate of side n is
     *   π(r,n) = (n mod (k_1 ··· k_r)) \ (k_1 ··· k_(r−1))
     * where k_i is the i^th component of the mesh dimensions.
     * See Rectangular_Meshes.pdf document for the derivation.
     */
    assert!(*r < *self.space_dim);
    assert!(*fe < self.num_fes);
    let cumprods_preceeding_ldims = if *r == 0 { 1 } else { self.cumprods_mesh_ldims[*r-1] };
    MeshCoord((*fe % self.cumprods_mesh_ldims[*r]) / cumprods_preceeding_ldims)
  }
  
  #[inline]
  fn fe_coord_mins_corner(&self, fe: FENum) -> ~[R] {
    vec::from_fn(*self.space_dim, |r| {
      let fe_mesh_coord_r = *self.fe_mesh_coord(Dim(r), fe) as R;
      self.min_bounds[r] + fe_mesh_coord_r * self.fe_dims[r]
    })
  }

}


impl<M:Monomial+RectIntegrable> Mesh<M>
                            for RectMesh<M> {

  #[inline(always)]
  fn num_fes(&self) -> uint {
    self.num_fes
  }

  #[inline(always)]
  fn num_nb_sides(&self) -> uint {
    self.num_nb_sides
  }
  
  #[inline(always)]
  fn num_oriented_element_shapes(&self) -> uint {
    1u
  }
  
  #[inline(always)]
  fn oriented_shape_for_fe(&self, fe: FENum) -> OShape {
    OShape(0)
  }
  
  #[inline(always)]
  fn num_side_faces_for_fe(&self, fe: FENum) -> uint {
    self.num_side_faces_per_fe
  }
  
  #[inline(always)]
  fn num_side_faces_for_shape(&self, oshape: OShape) -> uint {
    self.num_side_faces_per_fe 
  }
  
  #[inline(always)]
  fn dependent_dim_for_oshape_side(&self, oshape: OShape, side_face: SideFace) -> Dim {
    side_face_perp_axis(side_face)
  }
  
  #[inline]
  fn fe_inclusions_of_nb_side(&self, sn: NBSideNum) -> NBSideInclusions {
    let side_geom = self.nb_side_geom(sn);
    let a = side_geom.perp_axis;
    let lesser_fe = self.fe_with_mesh_coords(side_geom.mesh_coords);
    let greater_fe = FENum(*lesser_fe + (if *a == 0 {1} else {self.cumprods_mesh_ldims[*a-1]}));
    NBSideInclusions {
      nb_side_num: sn,
      fe1: lesser_fe,  sideface_in_fe1: greater_side_face_perp_to_axis(a),
      fe2: greater_fe, sideface_in_fe2: lesser_side_face_perp_to_axis(a) }
  }
  
  #[inline]
  fn nb_side_num_for_fe_side(&self, fe: FENum, side_face: SideFace) -> NBSideNum {
    let a = side_face_perp_axis(side_face);
    let side_mesh_coords = {
      if side_face_is_lesser_on_perp_axis(side_face) {
        let mut coords = self.fe_mesh_coords(fe);
        // lesser side perp axis coord = one less than corresponding fe coord
        coords[*a] = MeshCoord(*coords[*a]-1); 
        coords
      }
      else {
        self.fe_mesh_coords(fe)
      }
    };
    self.nb_side_with_mesh_coords(side_mesh_coords, a)
  }

  #[inline]
  fn is_boundary_side(&self, fe: FENum, side_face: SideFace) -> bool {
    let a = side_face_perp_axis(side_face);
    let coord_a = self.fe_mesh_coord(a, fe);
    let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
    *coord_a == 0 && is_lesser_side || !is_lesser_side && *coord_a == *self.mesh_ldims[*a]-1
  }
  
  fn num_boundary_sides(&self) -> uint {
    let d = *self.space_dim;
    let mut bsides = 0u;
    for perp_axis in range(0, d) {
      let mut prod = 1u;
      for r in range(0, d) {
        prod *= if r == perp_axis { 2 } else { *self.mesh_ldims[r] };
      }
      bsides += prod;
    }
    bsides
  }
  
  #[inline(always)]
  fn shape_diameter_inv(&self, oshape: OShape) -> R {
    self.rect_diameter_inv    
  }

  #[inline(always)]
  fn max_fe_diameter(&self) -> R {
    self.rect_diameter    
  }
  
  #[inline(always)]
  fn fe_interior_origin(&self, fe: FENum) -> ~[R] {
    self.fe_coord_mins_corner(fe)
  }
  
  #[inline]
  fn num_non_boundary_sides_for_fe(&self, fe: FENum) -> uint {
    range(0u8, self.num_side_faces_for_fe(fe) as u8)
      .count(|sf| self.is_boundary_side(fe, SideFace(sf))) 
  }
  
  #[inline]
  fn max_num_shape_sides(&self) -> uint {
    self.num_side_faces_per_fe 
  }


  // integration functions

  #[inline]
  fn intg_global_fn_on_fe_int(&self, f: &fn(&[R]) -> R, fe: FENum) -> R {
    let d = *self.space_dim;
    let fe_min_corner = self.fe_coord_mins_corner(fe);
    let fe_max_corner = vec::from_fn(d, |r| fe_min_corner[r] + self.fe_dims[r]);

    quadrature(&f, fe_min_corner, fe_max_corner, self.integration_rel_err, self.integration_abs_err)
  }

  #[inline]
  fn intg_global_fn_x_facerel_mon_on_fe_int(&self, f: &fn(&[R]) -> R, mon: M, fe: FENum) -> R {
    let d = *self.space_dim;
    let fe_min_corner = &self.fe_coord_mins_corner(fe);
    let fe_int_origin = fe_min_corner;
    let fe_max_corner = vec::from_fn(d, |r| fe_min_corner[r] + self.fe_dims[r]);
    
    quadrature(&|x: &[R]| { f(x) * mon.value_at_for_origin(x, *fe_int_origin) },
               *fe_min_corner, fe_max_corner,
               self.integration_rel_err, self.integration_abs_err)
  }

  #[inline]
  fn intg_facerel_poly_on_oshape_int<P:Polynomial<M>>(&self, p: P, oshape: OShape) -> R {
    p.foldl_terms(0 as R, |sum, (coef, mon)| {
      sum + coef * mon.integral_over_rect_at_origin(self.fe_dims)  
    })
  }

  #[inline]
  fn intg_facerel_poly_x_facerel_poly_on_oshape_int<P:Polynomial<M>>(&self, p1: P, p2: P, oshape: OShape) -> R {
    p1.foldl_terms(0 as R, |sum, (coef1, mon1)| {
      p2.foldl_terms(sum, |sum, (coef2, mon2)| {
        sum + coef1 * coef2 * (mon1*mon2).integral_over_rect_at_origin(self.fe_dims)
      })
    })
  }

  #[inline]
  fn intg_facerel_poly_x_facerel_poly_on_oshape_side<P:Polynomial<M>>(&self, p1: P, p2: P, oshape: OShape, side_face: SideFace) -> R {
    let a = side_face_perp_axis(side_face);
    p1.foldl_terms(0 as R, |sum, (coef1, mon1)| {
      p2.foldl_terms(sum, |sum, (coef2, mon2)| {
        sum + coef1 * coef2 * (mon1*mon2).surface_integral_siderel_over_rect_side(self.fe_dims, a)
      })
    })
  }

  #[inline]
  fn intg_facerel_mon_x_facerel_mon_on_oshape_int(&self, mon1: M, mon2: M, oshape: OShape) -> R {
    (mon1*mon2).integral_over_rect_at_origin(self.fe_dims)
  }

  #[inline]
  fn intg_facerel_mon_x_facerel_mon_on_oshape_side(&self, mon1: M, mon2: M, oshape: OShape, side_face: SideFace) -> R {
    let a = side_face_perp_axis(side_face);
    (mon1*mon2).surface_integral_siderel_over_rect_side(self.fe_dims, a)
  }

  #[inline]
  fn intg_facerel_mon_x_facerel_poly_on_oshape_int<P:Polynomial<M>>(&self, mon: M, p: P, oshape: OShape) -> R {
    p.foldl_terms(0 as R, |sum, (coef, p_mon)| {
      sum + coef * (mon*p_mon).integral_over_rect_at_origin(self.fe_dims)
    })
  }

  #[inline]
  fn intg_facerel_mon_x_facerel_poly_on_oshape_side<P:Polynomial<M>>(&self, mon: M, p: P, oshape: OShape, side_face: SideFace) -> R {
    let a = side_face_perp_axis(side_face);
    p.foldl_terms(0 as R, |sum, (coef, p_mon)| {
      sum + coef * (mon*p_mon).surface_integral_siderel_over_rect_side(self.fe_dims, a)
    })
  }

  #[inline]
  fn intg_intrel_mon_x_siderel_mon_on_oshape_side(&self, int_mon: M, side_mon: M, oshape: OShape, side_face: SideFace) -> R {
    let a = side_face_perp_axis(side_face);
    let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
    let side_intrel_a_coord = if is_lesser_side { 0 as R } else { self.fe_dims[*a] };

    /* Here we break the interior-relative monomial on the side into the constant a-dim factor and the monomial
       of other dimension factors. Since the interior and side-relative coordinate systems differ only in dimension 
       a, the latter monomial has the same expression in side-relative coordinates. */
    let int_mon_dim_a_fac = pow(side_intrel_a_coord, *int_mon.exp(a) as uint);
    let int_mon_wo_dim_a_fac = int_mon.map_exp(a, |_| Deg(0));

    int_mon_dim_a_fac * (int_mon_wo_dim_a_fac * side_mon).surface_integral_siderel_over_rect_side(self.fe_dims, a)
  }

  #[inline]
  fn intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(&self, side_mon: M, int_vmon: VectorMonomial<M>, oshape: OShape, side_face: SideFace) -> R {
    let a = side_face_perp_axis(side_face);
    match int_vmon.mon_dim() {
      Dim(r) if r == *a => {
        let int_vmon_mon = int_vmon.mon();
        let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
        let side_intrel_a_coord = if is_lesser_side { 0 as R } else { self.fe_dims[*a] };

        /* Here we break the interior-relative monomial on the side into the constant a-dim factor and the monomial
           of other dimension factors. Since the interior and side-relative coordinate systems differ only in dimension 
           a, the latter monomial has the same expression in side-relative coordinates. */
        let int_vmon_mon_dim_a_fac = pow(side_intrel_a_coord, *int_vmon_mon.exp(a) as uint);
        let int_vmon_mon_wo_dim_a_fac = int_vmon_mon.map_exp(a, |_| Deg(0));

        let outward_sense = if is_lesser_side { -1 as R } else { 1 as R };

        outward_sense * 
        int_vmon_mon_dim_a_fac *
        (int_vmon_mon_wo_dim_a_fac * side_mon).surface_integral_siderel_over_rect_side(self.fe_dims, a)
      }
      _ => 0 as R
    }
  }
 
  fn intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side<P:Polynomial<M>>(&self, p: P, int_vmon: VectorMonomial<M>, oshape: OShape, side_face: SideFace) -> R {
    let a = side_face_perp_axis(side_face);
    match int_vmon.mon_dim() {
      Dim(r) if r == *a => {
        let int_vmon_mon = int_vmon.mon();
        let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
        let side_intrel_a_coord = if is_lesser_side { 0 as R } else { self.fe_dims[*a] };

        /* Here we break the interior-relative monomial on the side into the constant a-dim factor and the monomial
           of other dimension factors. Since the interior and side-relative coordinate systems differ only in dimension 
           a, the latter monomial has the same expression in side-relative coordinates. */
        let int_vmon_mon_dim_a_fac = pow(side_intrel_a_coord, *int_vmon_mon.exp(a) as uint);
        let int_vmon_mon_wo_dim_a_fac = int_vmon_mon.map_exp(a, |_| Deg(0));

        let outward_sense = if is_lesser_side { -1 as R } else { 1 as R };

        outward_sense * 
        int_vmon_mon_dim_a_fac *
        p.foldl_terms(0 as R, |sum, (coef, mon)| {
          sum + coef * (int_vmon_mon_wo_dim_a_fac * mon).surface_integral_siderel_over_rect_side(self.fe_dims, a)
        })
      }
      _ => 0 as R
    }
  }

}

// side-related auxiliary stateless functions

// Find the axis which is perpendicular to the given side face.
#[inline]
fn side_face_perp_axis(side_face: SideFace) -> Dim {
  Dim(*side_face as uint / 2)
}

// Determine whether a side face is the one with lesser axis value along its perpendicular axis.
#[inline(always)]
fn side_face_is_lesser_on_perp_axis(side_face: SideFace) -> bool {
  *side_face % 2 == 0
}

// Returns the side face of lesser coordinate value along the indicated axis.
#[inline]
fn lesser_side_face_perp_to_axis(a: Dim) -> SideFace {
  SideFace((2 * *a) as u8)
}

// Returns the side face of greater coordinate value along the indicated axis.
#[inline]
fn greater_side_face_perp_to_axis(a: Dim) -> SideFace {
  SideFace((2 * *a + 1) as u8)
}


// integrable trait

trait RectIntegrable {

  /// Integrate the monomial over the rectangle of indicated dimensions having its minimums corner at the origin.
  fn integral_over_rect_at_origin(&self, rect_dims: &[R]) -> R;  

  /// Integrate a *side-relative* monomial over a side face of a rectangle of given dimenions, the side being
  /// perpendicular to the indicated axis. The monomial's origin is the side's corner of minimum coordinates.
  /// Note that the integral will have the same value no matter which of the two sides perpendicular to the
  /// indicated axis is chosen as the domain of integration, because of the side local coordinates used to
  /// interpet the monomial (the  monomial will either have 0 exponent for the perpendicular axis variable
  /// and thus be independent of values in that coordinate, or else the integral will be 0).
  fn surface_integral_siderel_over_rect_side(&self, rect_dims: &[R], side_perp_axis: Dim) -> R;

}

impl RectIntegrable for Mon1d {
  #[inline]
  fn integral_over_rect_at_origin(&self, rect_dims: &[R]) -> R {
    let exp_plus_1 = *self.exps[0] as uint + 1;
    pow(rect_dims[0], exp_plus_1)/(exp_plus_1 as R)
  }

  fn surface_integral_siderel_over_rect_side(&self, rect_dims: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[0] != 0u8 {0 as R} else {1 as R}
  }
}

impl RectIntegrable for Mon2d {
  #[inline]
  fn integral_over_rect_at_origin(&self, rect_dims: &[R]) -> R {
    let exp0_plus_1 = *self.exps[0] as uint + 1;
    let exp1_plus_1 = *self.exps[1] as uint + 1;
    pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
    pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R)
  }
  #[inline]
  fn surface_integral_siderel_over_rect_side(&self, rect_dims: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[*side_perp_axis] != 0u8 { 0 as R }
    else {
      match side_perp_axis {
        Dim(0) => {
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R)
        }
        Dim(1) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R)
        }
        _ => fail!("Dimension out of range.")
      }
    }
  }
}

impl RectIntegrable for Mon3d {
  #[inline]
  fn integral_over_rect_at_origin(&self, rect_dims: &[R]) -> R {
    let exp0_plus_1 = *self.exps[0] as uint + 1;
    let exp1_plus_1 = *self.exps[1] as uint + 1;
    let exp2_plus_1 = *self.exps[2] as uint + 1;
    pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
    pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R) *
    pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R)
  }
  #[inline]
  fn surface_integral_siderel_over_rect_side(&self, rect_dims: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[*side_perp_axis] != 0u8 { 0 as R }
    else {
      match side_perp_axis {
        Dim(0) => {
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R)
        }
        Dim(1) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R)
        }
        Dim(2) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R)
        }
        _ => fail!("Dimension out of range.")
      }
    }
  }
}

impl RectIntegrable for Mon4d {
  #[inline]
  fn integral_over_rect_at_origin(&self, rect_dims: &[R]) -> R {
    let exp0_plus_1 = *self.exps[0] as uint + 1;
    let exp1_plus_1 = *self.exps[1] as uint + 1;
    let exp2_plus_1 = *self.exps[2] as uint + 1;
    let exp3_plus_1 = *self.exps[3] as uint + 1;
    pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
    pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R) *
    pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R) *
    pow(rect_dims[3], exp3_plus_1)/(exp3_plus_1 as R)
  }
  #[inline]
  fn surface_integral_siderel_over_rect_side(&self, rect_dims: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[*side_perp_axis] != 0u8 { 0 as R }
    else {
      match side_perp_axis {
        Dim(0) => {
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          let exp3_plus_1 = *self.exps[3] as uint + 1;
          pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R) *
          pow(rect_dims[3], exp3_plus_1)/(exp3_plus_1 as R)
        }
        Dim(1) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          let exp3_plus_1 = *self.exps[3] as uint + 1;
          pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R) *
          pow(rect_dims[3], exp3_plus_1)/(exp3_plus_1 as R)
        }
        Dim(2) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp3_plus_1 = *self.exps[3] as uint + 1;
          pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(rect_dims[3], exp3_plus_1)/(exp3_plus_1 as R)
        }
        Dim(3) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          pow(rect_dims[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(rect_dims[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(rect_dims[2], exp2_plus_1)/(exp2_plus_1 as R)
        }
        _ => fail!("Dimension out of range.")
      }
    }
  }
}







// Tests




#[test]
fn test_3x4x5_constr() -> () {
  let mesh_min_coords = ~[1.0f64, 2.0f64, 3.0f64];
  let mesh_max_coords = ~[2.0f64, 3., 4.];
  let mesh_ldims = ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)];
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(mesh_min_coords.clone(), mesh_max_coords.clone(), mesh_ldims.clone());
  let rect_oshape = OShape(0);

  assert_eq!(rmesh3x4x5.space_dim, Dim(3));
  assert_eq!(&rmesh3x4x5.min_bounds, &mesh_min_coords);
  assert_eq!(&rmesh3x4x5.max_bounds, &mesh_max_coords);
  assert_eq!(&rmesh3x4x5.mesh_ldims, &mesh_ldims);
  assert_eq!(&rmesh3x4x5.fe_dims, &~[1./3., 1./4., 1./5.]);
  assert_approx(rmesh3x4x5.rect_diameter, sqrt(pow(1./3.,2) + pow(1./4.,2) + pow(1./5.,2)));
  assert_approx(rmesh3x4x5.shape_diameter_inv(rect_oshape), 1./sqrt(pow(1./3.,2) + pow(1./4.,2) + pow(1./5.,2)));
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
fn test_3x4x5_mesh_coords() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
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
#[should_fail]
fn test_3x4x5_bad_mesh_coords() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let bad_access = rmesh3x4x5.fe_mesh_coord(Dim(0), FENum(5*12));
}

fn test_3x4x5_boundary_side_determ() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_face = lesser_side_face_perp_to_axis(Dim(0));
  let right_face = greater_side_face_perp_to_axis(Dim(0));
  let bottom_face = lesser_side_face_perp_to_axis(Dim(1));
  let top_face = greater_side_face_perp_to_axis(Dim(1));
  let back_face = lesser_side_face_perp_to_axis(Dim(2));
  let front_face = greater_side_face_perp_to_axis(Dim(2));

  assert!( rmesh3x4x5.is_boundary_side(FENum(0), left_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(0), right_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(0), bottom_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(0), top_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(0), back_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(0), front_face));

  assert!(!rmesh3x4x5.is_boundary_side(FENum(2), left_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(2), right_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(2), bottom_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(2), top_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(2), back_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(2), front_face));

  assert!( rmesh3x4x5.is_boundary_side(FENum(3), left_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), right_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), bottom_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), top_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(3), back_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(3), front_face));

  assert!(!rmesh3x4x5.is_boundary_side(FENum(11), left_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(11), right_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(11), bottom_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(11), top_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(11), back_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(11), front_face));

  assert!( rmesh3x4x5.is_boundary_side(FENum(12), left_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), right_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(12), bottom_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), top_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), back_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(12), front_face));

  assert!( rmesh3x4x5.is_boundary_side(FENum(4*12), left_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(4*12), right_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(4*12), bottom_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(4*12), top_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(4*12), back_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(4*12), front_face));

  assert!(!rmesh3x4x5.is_boundary_side(FENum(5*12-1), left_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(5*12-1), right_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(5*12-1), bottom_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(5*12-1), top_face));
  assert!(!rmesh3x4x5.is_boundary_side(FENum(5*12-1), back_face));
  assert!( rmesh3x4x5.is_boundary_side(FENum(5*12-1), front_face));
}

#[test]
#[should_fail]
fn test_3x4x5_bad_is_boundary_side_fenum() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let front_face = greater_side_face_perp_to_axis(Dim(2));
  let bad = rmesh3x4x5.is_boundary_side(FENum(5*12), front_face);
}

// Test the non-boundary sides perpendicular to axis 0.
#[test]
fn test_3x4x5_nonboundary_side_coords_axis0() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_face = lesser_side_face_perp_to_axis(Dim(0));
  let right_face = greater_side_face_perp_to_axis(Dim(0));

  // The mesh for non-boundary sides perpendicular to axis 0 has dimensions 2 x 4 x 5.

  // first side in first row (axis 0)
  let sgeom_0 = rmesh3x4x5.nb_side_geom(NBSideNum(0));
  assert_eq!(sgeom_0.perp_axis, Dim(0));
  assert_eq!(&sgeom_0.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(0));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(0)),
             NBSideInclusions { nb_side_num: NBSideNum(0),
                                fe1: FENum(0), sideface_in_fe1: right_face,
                                fe2: FENum(1), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(0), right_face), NBSideNum(0));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(1), left_face),  NBSideNum(0));

  // second side in first row
  let sgeom_1 = rmesh3x4x5.nb_side_geom(NBSideNum(1));
  assert_eq!(sgeom_1.perp_axis, Dim(0));
  assert_eq!(&sgeom_1.mesh_coords, &~[MeshCoord(1), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(0), MeshCoord(0)], Dim(0)), NBSideNum(1));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(1)),
             NBSideInclusions { nb_side_num: NBSideNum(1),
                                fe1: FENum(1), sideface_in_fe1: right_face,
                                fe2: FENum(2), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(1), right_face), NBSideNum(1));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(2), left_face),  NBSideNum(1));

  // first side in second row
  let sgeom_2 = rmesh3x4x5.nb_side_geom(NBSideNum(2));
  assert_eq!(sgeom_2.perp_axis, Dim(0));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(0)), NBSideNum(2));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(2)),
             NBSideInclusions { nb_side_num: NBSideNum(2),
                                fe1: FENum(3), sideface_in_fe1: right_face,
                                fe2: FENum(4), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3), right_face), NBSideNum(2));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(4), left_face),  NBSideNum(2));

  // last side in first stack
  let sgeom_7 = rmesh3x4x5.nb_side_geom(NBSideNum(7));
  assert_eq!(sgeom_7.perp_axis, Dim(0));
  assert_eq!(&sgeom_7.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(0)], Dim(0)), NBSideNum(7));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(7)),
             NBSideInclusions { nb_side_num: NBSideNum(7),
                                fe1: FENum(10), sideface_in_fe1: right_face,
                                fe2: FENum(11), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(10), right_face), NBSideNum(7));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(11), left_face),  NBSideNum(7));

  // first side of second stack
  let sgeom_8 = rmesh3x4x5.nb_side_geom(NBSideNum(8));
  assert_eq!(sgeom_8.perp_axis, Dim(0));
  assert_eq!(&sgeom_8.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(0)), NBSideNum(8));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(8)),
             NBSideInclusions { nb_side_num: NBSideNum(8),
                                fe1: FENum(12), sideface_in_fe1: right_face,
                                fe2: FENum(13), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), right_face), NBSideNum(8));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(13), left_face),  NBSideNum(8));

  // first side of last stack
  let sgeom_32 = rmesh3x4x5.nb_side_geom(NBSideNum(32));
  assert_eq!(sgeom_32.perp_axis, Dim(0));
  assert_eq!(&sgeom_32.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4)], Dim(0)), NBSideNum(32));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(32)),
             NBSideInclusions { nb_side_num: NBSideNum(32),
                                fe1: FENum(48), sideface_in_fe1: right_face,
                                fe2: FENum(49), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), right_face), NBSideNum(32));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(49), left_face),  NBSideNum(32));

  // last side of last stack
  let sgeom_39 = rmesh3x4x5.nb_side_geom(NBSideNum(39));
  assert_eq!(sgeom_39.perp_axis, Dim(0));
  assert_eq!(&sgeom_39.mesh_coords, &~[MeshCoord(1), MeshCoord(3), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(1), MeshCoord(3), MeshCoord(4)], Dim(0)), NBSideNum(39));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(39)),
             NBSideInclusions { nb_side_num: NBSideNum(39),
                                fe1: FENum(58), sideface_in_fe1: right_face,
                                fe2: FENum(59), sideface_in_fe2: left_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(58), right_face), NBSideNum(39));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(59), left_face),  NBSideNum(39));
}

// Test the non-boundary sides perpendicular to axis 1.
#[test]
fn test_3x4x5_nonboundary_side_coords_axis1() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let bottom_face = lesser_side_face_perp_to_axis(Dim(1));
  let top_face = greater_side_face_perp_to_axis(Dim(1));

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
                                fe1: FENum(0), sideface_in_fe1: top_face,
                                fe2: FENum(3), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(0), top_face), NBSideNum(first_axis1+0));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3), bottom_face),  NBSideNum(first_axis1+0));

  // last side in first row
  let sgeom_2 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+2));
  assert_eq!(sgeom_2.perp_axis, Dim(1));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(2), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(0), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+2));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+2)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+2),
                                fe1: FENum(2), sideface_in_fe1: top_face,
                                fe2: FENum(5), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(2), top_face), NBSideNum(first_axis1+2));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(5), bottom_face),  NBSideNum(first_axis1+2));

  // first side in second row
  let sgeom_3 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+3));
  assert_eq!(sgeom_3.perp_axis, Dim(1));
  assert_eq!(&sgeom_3.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+3));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+3)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+3),
                                fe1: FENum(3), sideface_in_fe1: top_face,
                                fe2: FENum(6), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3), top_face), NBSideNum(first_axis1+3));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(6), bottom_face),  NBSideNum(first_axis1+3));

  // last side in first stack
  let sgeom_8 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+8));
  assert_eq!(sgeom_8.perp_axis, Dim(1));
  assert_eq!(&sgeom_8.mesh_coords, &~[MeshCoord(2), MeshCoord(2), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(2), MeshCoord(0)], Dim(1)), NBSideNum(first_axis1+8));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+8)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+8),
                                fe1: FENum(8), sideface_in_fe1: top_face,
                                fe2: FENum(11), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(8), top_face), NBSideNum(first_axis1+8));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(11), bottom_face),  NBSideNum(first_axis1+8));

  // first side of second stack
  let sgeom_9 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+9));
  assert_eq!(sgeom_9.perp_axis, Dim(1));
  assert_eq!(&sgeom_9.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(1)), NBSideNum(first_axis1+9));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+9)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+9),
                                fe1: FENum(12), sideface_in_fe1: top_face,
                                fe2: FENum(15), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), top_face), NBSideNum(first_axis1+9));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(15), bottom_face),  NBSideNum(first_axis1+9));

  // first side of last stack
  let sgeom_36 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+36));
  assert_eq!(sgeom_36.perp_axis, Dim(1));
  assert_eq!(&sgeom_36.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4)], Dim(1)), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+36)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+36),
                                fe1: FENum(48), sideface_in_fe1: top_face,
                                fe2: FENum(51), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), top_face), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(51), bottom_face),  NBSideNum(first_axis1+36));

  // last side of last stack
  let sgeom_36 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis1+36));
  assert_eq!(sgeom_36.perp_axis, Dim(1));
  assert_eq!(&sgeom_36.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(4)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(4)], Dim(1)), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis1+36)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis1+36),
                                fe1: FENum(48), sideface_in_fe1: top_face,
                                fe2: FENum(51), sideface_in_fe2: bottom_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), top_face), NBSideNum(first_axis1+36));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(51), bottom_face),  NBSideNum(first_axis1+36));

}

// Test the non-boundary sides perpendicular to axis 2.
#[test]
fn test_3x4x5_nonboundary_side_coords_axis2() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let front_face = lesser_side_face_perp_to_axis(Dim(2));
  let back_face = greater_side_face_perp_to_axis(Dim(2));

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
                                fe1: FENum(0),  sideface_in_fe1: back_face,
                                fe2: FENum(12), sideface_in_fe2: front_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(0),  back_face), NBSideNum(first_axis2+0));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), front_face),  NBSideNum(first_axis2+0));

  // third side in first row
  let sgeom_2 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+2));
  assert_eq!(sgeom_2.perp_axis, Dim(2));
  assert_eq!(&sgeom_2.mesh_coords, &~[MeshCoord(2), MeshCoord(0), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(0), MeshCoord(0)], Dim(2)), NBSideNum(first_axis2+2));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+2)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+2),
                                fe1: FENum(2),  sideface_in_fe1: back_face,
                                fe2: FENum(14), sideface_in_fe2: front_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(2),  back_face), NBSideNum(first_axis2+2));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(14), front_face),  NBSideNum(first_axis2+2));


  // second row
  let sgeom_3 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+3));
  assert_eq!(sgeom_3.perp_axis, Dim(2));
  assert_eq!(&sgeom_3.mesh_coords, &~[MeshCoord(0), MeshCoord(1), MeshCoord(0)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(1), MeshCoord(0)], Dim(2)), NBSideNum(first_axis2+3));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+3)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+3),
                                fe1: FENum(3),  sideface_in_fe1: back_face,
                                fe2: FENum(15), sideface_in_fe2: front_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(3),  back_face), NBSideNum(first_axis2+3));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(15), front_face),  NBSideNum(first_axis2+3));
  
  // first side of second stack
  let sgeom_12 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+12));
  assert_eq!(sgeom_12.perp_axis, Dim(2));
  assert_eq!(&sgeom_12.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(1)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(1)], Dim(2)), NBSideNum(first_axis2+12));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+12)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+12),
                                fe1: FENum(12), sideface_in_fe1: back_face,
                                fe2: FENum(24), sideface_in_fe2: front_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(12), back_face), NBSideNum(first_axis2+12));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(24), front_face),  NBSideNum(first_axis2+12));

  // first side of last stack
  let sgeom_36 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+36));
  assert_eq!(sgeom_36.perp_axis, Dim(2));
  assert_eq!(&sgeom_36.mesh_coords, &~[MeshCoord(0), MeshCoord(0), MeshCoord(3)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(0), MeshCoord(0), MeshCoord(3)], Dim(2)), NBSideNum(first_axis2+36));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+36)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+36),
                                fe1: FENum(36), sideface_in_fe1: back_face,
                                fe2: FENum(48), sideface_in_fe2: front_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(36), back_face), NBSideNum(first_axis2+36));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(48), front_face),  NBSideNum(first_axis2+36));

  // last side of last stack
  let sgeom_47 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+47));
  assert_eq!(sgeom_47.perp_axis, Dim(2));
  assert_eq!(&sgeom_47.mesh_coords, &~[MeshCoord(2), MeshCoord(3), MeshCoord(3)]);
  assert_eq!(rmesh3x4x5.nb_side_with_mesh_coords(&[MeshCoord(2), MeshCoord(3), MeshCoord(3)], Dim(2)), NBSideNum(first_axis2+47));
  assert_eq!(rmesh3x4x5.fe_inclusions_of_nb_side(NBSideNum(first_axis2+47)),
             NBSideInclusions { nb_side_num: NBSideNum(first_axis2+47),
                                fe1: FENum(47), sideface_in_fe1: back_face,
                                fe2: FENum(59), sideface_in_fe2: front_face });
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(47), back_face), NBSideNum(first_axis2+47));
  assert_eq!(rmesh3x4x5.nb_side_num_for_fe_side(FENum(59), front_face), NBSideNum(first_axis2+47));

  // Above was the last non-boundary side.  Check that the total number is correct.
  assert_eq!(rmesh3x4x5.num_boundary_sides(), 2*20 + 2*15 + 2*12);
}

// Attempt to access beyond the last non-boundary side should fail.
#[test]
#[should_fail]
fn test_3x4x5_nonboundary_bad_side_coords_axis2() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let first_axis2 = 2*4*5 + 3*3*5;
  let sgeom_47 = rmesh3x4x5.nb_side_geom(NBSideNum(first_axis2+48));
}

#[test]
fn test_3x4x5_conv_fe_coords_to_fenum() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(0)]), FENum(0));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(2), MeshCoord(0), MeshCoord(0)]), FENum(2));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(0), MeshCoord(1), MeshCoord(0)]), FENum(3));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(1), MeshCoord(1), MeshCoord(0)]), FENum(4));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(0), MeshCoord(0), MeshCoord(1)]), FENum(12));
  assert_eq!(rmesh3x4x5.fe_with_mesh_coords([MeshCoord(2), MeshCoord(3), MeshCoord(4)]), FENum(59));
}

// Test integration of monomials through RectIntegrable trait.

#[test]
fn test_intg_1d() {
  let one = Mon1d { exps: [Deg(0)] };
  let x = Mon1d { exps: [Deg(1)] };
  let x2 = Mon1d { exps: [Deg(2)] };
  assert_eq!(one.integral_over_rect_at_origin([2.]), 2.);
  assert_eq!(x.integral_over_rect_at_origin([2.]), 2.);
  assert_eq!(x2.integral_over_rect_at_origin([2.]), 8./3.);
}

#[test]
fn test_intg_2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x1y2 = Mon2d { exps: [Deg(1), Deg(2)] };
  let x3y1 = Mon2d { exps: [Deg(3), Deg(1)] };
  assert_eq!(one.integral_over_rect_at_origin([2.,3.]), 6.);
  assert_eq!(x1y2.integral_over_rect_at_origin([2.,3.]), 18.);
  assert_eq!(x3y1.integral_over_rect_at_origin([2.,3.]), 4.*9./2.);
}

#[test]
fn test_intg_3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x1y2z3 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  let x3y1z2 = Mon3d { exps: [Deg(3), Deg(1), Deg(2)] };
  assert_eq!(one.integral_over_rect_at_origin([2.,3.,4.]), 24.);
  assert_eq!(x1y2z3.integral_over_rect_at_origin([2.,3.,4.]), 18.*64.);
  assert_eq!(x3y1z2.integral_over_rect_at_origin([2.,3.,4.]), 4.*(9./2.)*64./3.);
}

#[test]
fn test_intg_4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  let x1y2z3t4 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  let x3y1z2t4 = Mon4d { exps: [Deg(3), Deg(1), Deg(2), Deg(4)] };
  assert_eq!(one.integral_over_rect_at_origin([2.,3.,4.,5.]), 24.*5.);
  assert_eq!(x1y2z3t4.integral_over_rect_at_origin([2.,3.,4.,5.]), 18.*64.*625.);
  assert_eq!(x3y1z2t4.integral_over_rect_at_origin([2.,3.,4.,5.]), 4.*(9./2.)*(64./3.)*625.);
}


#[test]
fn test_intg_const_facerel_mons() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  let left_face = lesser_side_face_perp_to_axis(Dim(0));
  let right_face = greater_side_face_perp_to_axis(Dim(0));
  let bottom_face = lesser_side_face_perp_to_axis(Dim(1));
  let top_face = greater_side_face_perp_to_axis(Dim(1));
  let back_face = lesser_side_face_perp_to_axis(Dim(2));
  let front_face = greater_side_face_perp_to_axis(Dim(2));

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_int(one, one, OShape(0)),
                1./3. * 1./4. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(one, one, OShape(0), left_face),
                1./4. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(one, one, OShape(0), right_face),
                1./4. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(one, one, OShape(0), bottom_face),
                1./3. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(one, one, OShape(0), top_face),
                1./3. * 1./5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(one, one, OShape(0), back_face),
                1./3. * 1./4.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(one, one, OShape(0), front_face),
                1./3. * 1./4.);
}

#[test]
fn test_intg_global_fn_on_fe_int() -> () {
  let mins = ~[1.0f64, 2.0f64, 3.0f64];
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(mins.clone(),
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x3_y4_z = |x:&[R]| -> R {
    pow(x[0]-mins[0], 3) * pow(x[1]-mins[1], 4) * (x[2]-mins[2])
  };
  
  assert_approx(rmesh3x4x5.intg_global_fn_on_fe_int(x3_y4_z, FENum(0)),
                pow(1./3.,4)/4. * pow(1./4.,5)/5. * pow(1./5.,2)/2.);
}

#[test]
fn test_intg_global_fn_x_mon_on_fe0_int() -> () {
  let mins = ~[1.0f64, 2.0f64, 3.0f64];
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(mins.clone(),
                                                   ~[2.0f64, 3., 4.],
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
fn test_intg_global_fn_x_mon_on_fe16_int() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
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
fn test_intg_mons_on_fe_int() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_int(x*y*z, z, OShape(0)),
                pow(1./3.,2)/2. * pow(1./4.,2)/2. * pow(1./5.,3)/3.);
}

#[test]
fn test_intg_global_fn_x_mon_on_fe_sides() -> () {
  let rmesh3x4x5: ~RectMesh<Mon3d> = RectMesh::new(~[1.0f64, 2.0f64, 3.0f64],
                                                   ~[2.0f64, 3., 4.],
                                                   ~[MeshCoord(3), MeshCoord(4), MeshCoord(5)]);
  let left_face = lesser_side_face_perp_to_axis(Dim(0));
  let right_face = greater_side_face_perp_to_axis(Dim(0));
  let bottom_face = lesser_side_face_perp_to_axis(Dim(1));
  let top_face = greater_side_face_perp_to_axis(Dim(1));
  let back_face = lesser_side_face_perp_to_axis(Dim(2));
  let front_face = greater_side_face_perp_to_axis(Dim(2));

  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(y*z, z, OShape(0), left_face),
                pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(y*z, z, OShape(0), right_face),
                pow(1./4.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(y*z, z, OShape(0), bottom_face),
                0.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(y*z, z, OShape(0), top_face),
                0.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*z, z, OShape(0), bottom_face),
                pow(1./3.,2)/2. * pow(1./5.,3)/3.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*x, y, OShape(0), back_face),
                pow(1./3.,3)/3. * pow(1./4.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*x, y, OShape(0), front_face),
                pow(1./3.,3)/3. * pow(1./4.,2)/2.);

  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(y*y, y*y*z, OShape(0), left_face),
                pow(1./4.,5)/5. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(y*y, y*y*z, OShape(0), right_face),
                pow(1./4.,5)/5. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*x, x*z, OShape(0), bottom_face),
                pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*x, x*z, OShape(0), top_face),
                pow(1./3.,4)/4. * pow(1./5.,2)/2.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*x, x*y*y*y*y, OShape(0), back_face),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
  assert_approx(rmesh3x4x5.intg_facerel_mon_x_facerel_mon_on_oshape_side(x*x, x*y*y*y*y, OShape(0), front_face),
                pow(1./3.,4)/4. * pow(1./4.,5)/5.);
}


/*
# Test integrals on finite element faces.

# Integrate monomials on faces of reference element.

@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(x^3*y^4*z, rect_oshape, Mesh.interior_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2
)


# Integrate polynomial on faces of reference element.
@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(2*y*z^2 + 3*x^3*y^4*z, rect_oshape, Mesh.interior_face, rmesh3x4x5),
  2(1/3 * (1/4)^2/2 * (1/5)^3/3) + 3((1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(1.2*y*z^2 + 3.4*x^3*y^4*z, rect_oshape, left_face, rmesh3x4x5),
  1.2((1/4)^2/2 * (1/5)^3/3)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(1.5*y*z^2 + 0.123*y^4*z, rect_oshape, right_face, rmesh3x4x5),
  1.5((1/4)^2/2 * (1/5)^3/3) + 0.123((1/4)^5/5 * (1/5)^2/2)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(y*z^2 + x^3*y^4*z, rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(4.5*x*z^2 + 2.3x^3*z, rect_oshape, top_face, rmesh3x4x5),
  4.5((1/3)^2/2 * (1/5)^3/3) + 2.3((1/3)^4/4 * (1/5)^2/2)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_oshape_face(x*y^2 + x^3*y^4, rect_oshape, back_face, rmesh3x4x5),
  (1/3)^2/2 * (1/4)^3/3 + (1/3)^4/4 * (1/4)^5/5
)


# Integrate vector monomials vs outward normals on the side faces of the reference finite element.

@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, right_face, rmesh3x4x5),
  (1/3)^3 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(z, VM(y^4*z, dim(1)), rect_oshape, left_face, rmesh3x4x5),
  -(1/4)^5/5 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y^2*z,   VM(x^3*y^4*z, dim(1)), rect_oshape, right_face, rmesh3x4x5),
  (1/3)^3 * (1/4)^7/7 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x*y^2*z, VM(x^3*y^4*z, dim(1)), rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, front_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, top_face, rmesh3x4x5),
  (1/4)^4 * (1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x^2*z,   VM(x^3*y^4*z, dim(2)), rect_oshape, top_face, rmesh3x4x5),
  (1/4)^4 * (1/3)^6/6 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x^2*y*z, VM(x^3*y^4*z, dim(2)), rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, front_face, rmesh3x4x5),
  0
)

@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^4/4 * (1/4)^5/5
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x^2*y,   VM(x^3*y^4*z, dim(3)), rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^6/6 * (1/4)^6/6
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x^2*y*z, VM(x^3*y^4*z, dim(3)), rect_oshape, front_face, rmesh3x4x5),
  0
)

# Integrate on lesser faces using monomials which are constant in the corresponding dimension (so the integrals aren't 0).
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(y^4*z,   dim(1)), rect_oshape, left_face, rmesh3x4x5),
  -(1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*z,   dim(2)), rect_oshape, bottom_face, rmesh3x4x5),
  -(1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4, dim(3)), rect_oshape, back_face, rmesh3x4x5),
  -(1/3)^4/4 * (1/4)^5/5
)


# global vs polynomial
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, 2x*y*z + 4.5x*y*z, fenum(17), Mesh.interior_face, rmesh3x4x5),
  (2 + 4.5)*(1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2
)


@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*y, y, rect_oshape, right_face, rmesh3x4x5),
  (1/3) * (1/4)^3/3 * (1/5)
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*y*z, y*z, rect_oshape, right_face, rmesh3x4x5),
  (1/3) * (1/4)^3/3 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y*z, x*y*z, rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*y, y, rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y^2*z, y*z, rect_oshape, left_face, rmesh3x4x5),
  (1/4)^4/4 * (1/5)^3/3
)

@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y*z, y, rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y*z, x, rect_oshape, top_face, rmesh3x4x5),
  (1/4) * (1/3)^2/2 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*y*z, x*z, rect_oshape, top_face, rmesh3x4x5),
  (1/4) * (1/3)^3/3 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*y, x*z, rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*z, z, rect_oshape, bottom_face, rmesh3x4x5),
  (1/3)^2/2 * (1/5)^3/3
)

@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y*x, z, rect_oshape, front_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(z*x, y, rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^2/2 * (1/4)^2/2
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(z*x*y, x*y, rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^3/3 * (1/4)^3/3
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y*z, y*x, rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x*y, x*y^2, rect_oshape, back_face, rmesh3x4x5),
  (1/3)^3/3 * (1/4)^4/4
)

# 10 rows x 20 cols mesh, with vertexes at each integer pair

rmesh20x10 = RectMesh([0.,0.], [20.,10.], [MeshCoord(20), MeshCoord(10)])

@test Mesh.fe_interior_origin(fenum(1), rmesh20x10) == [0., 0.]
@test Mesh.fe_interior_origin(fenum(2), rmesh20x10) == [1., 0.]
@test Mesh.fe_interior_origin(fenum(20), rmesh20x10) == [19., 0.]
@test Mesh.fe_interior_origin(fenum(21), rmesh20x10) == [0., 1.]
@test Mesh.fe_interior_origin(fenum(200), rmesh20x10) == [19., 9.]

@test Mesh.num_fes(rmesh20x10) == 200
@test Mesh.num_nb_sides(rmesh20x10) == 370

# test boundary sides

@test Mesh.num_boundary_sides(rmesh20x10) == 2*10 + 2*20

@test !Mesh.is_boundary_side(fenum(1), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(1), right_face, rmesh20x10)
@test  Mesh.is_boundary_side(fenum(1), bottom_face,  rmesh20x10)
@test  Mesh.is_boundary_side(fenum(1), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fenum(8), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(8), right_face, rmesh20x10)
@test  Mesh.is_boundary_side(fenum(8), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fenum(8), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fenum(20), top_face, rmesh20x10)
@test  Mesh.is_boundary_side(fenum(20), right_face, rmesh20x10)
@test  Mesh.is_boundary_side(fenum(20), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fenum(20), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fenum(21), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(21), right_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(21), bottom_face,  rmesh20x10)
@test  Mesh.is_boundary_side(fenum(21), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fenum(22), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(22), right_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(22), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fenum(22), left_face,  rmesh20x10)

@test  Mesh.is_boundary_side(fenum(200), top_face, rmesh20x10)
@test  Mesh.is_boundary_side(fenum(200), right_face, rmesh20x10)
@test !Mesh.is_boundary_side(fenum(200), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fenum(200), left_face,  rmesh20x10)


rmesh3x2 = RectMesh([0.,0.], [3.,2.], [MeshCoord(3), MeshCoord(2)])

@test Mesh.num_boundary_sides(rmesh3x2) == 2*3 + 2*2

@test RMesh.perp_axis_for_nb_side(nbsidenum(1), rmesh3x2) == dim(1)
@test RMesh.perp_axis_for_nb_side(nbsidenum(4), rmesh3x2) == dim(1)
@test RMesh.perp_axis_for_nb_side(nbsidenum(5), rmesh3x2) == dim(2)
@test RMesh.perp_axis_for_nb_side(nbsidenum(7), rmesh3x2) == dim(2)
@test_throws RMesh.perp_axis_for_nb_side(nbsidenum(8), rmesh3x2)

# Test side inclusions
# fe vertical sides
incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(1), rmesh3x2)
@test incls.fe1 == fenum(1)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fenum(2)
@test incls.face_in_fe2 == left_face

incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(2), rmesh3x2)
@test incls.fe1 == fenum(2)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fenum(3)
@test incls.face_in_fe2 == left_face

incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(3), rmesh3x2)
@test incls.fe1 == fenum(4)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fenum(5)
@test incls.face_in_fe2 == left_face

incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(4), rmesh3x2)
@test incls.fe1 == fenum(5)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fenum(6)
@test incls.face_in_fe2 == left_face

# fe horizontal sides
incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(5), rmesh3x2)
@test incls.fe1 == fenum(1)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fenum(4)
@test incls.face_in_fe2 == bottom_face

incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(6), rmesh3x2)
@test incls.fe1 == fenum(2)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fenum(5)
@test incls.face_in_fe2 == bottom_face

incls = Mesh.fe_inclusions_of_nb_side(nbsidenum(7), rmesh3x2)
@test incls.fe1 == fenum(3)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fenum(6)
@test incls.face_in_fe2 == bottom_face


# Test integrals on finite element faces.

x = Monomial(1,0)
y = Monomial(0,1)
one_mon = RMesh.one_mon(rmesh3x2)

# Integrate constants on faces of reference element.
@test Mesh.integral_face_rel_on_oshape_face(one_mon, rect_oshape, Mesh.interior_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_oshape_face(one_mon, rect_oshape, left_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_oshape_face(one_mon, rect_oshape, right_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_oshape_face(one_mon, rect_oshape, bottom_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_oshape_face(one_mon, rect_oshape, top_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_oshape_face(3., rect_oshape, Mesh.interior_face, rmesh3x2) == 3.
@test Mesh.integral_face_rel_on_oshape_face(3., rect_oshape, top_face, rmesh3x2) == 3.

# Integrate monomials on faces of reference element.

@test Mesh.integral_face_rel_on_oshape_face(x*y^2, rect_oshape, Mesh.interior_face, rmesh3x2) == 1/6
@test Mesh.integral_face_rel_on_oshape_face(y^2, rect_oshape, left_face, rmesh3x2) == 1/3
@test Mesh.integral_face_rel_on_oshape_face(y^2, rect_oshape, right_face, rmesh3x2) == 1/3
@test Mesh.integral_face_rel_on_oshape_face(x*y^2, rect_oshape, right_face, rmesh3x2) == 0
@test Mesh.integral_face_rel_on_oshape_face(x, rect_oshape, bottom_face, rmesh3x2) == 1/2
@test Mesh.integral_face_rel_on_oshape_face(x, rect_oshape, top_face, rmesh3x2) == 1/2
@test Mesh.integral_face_rel_on_oshape_face(x*y, rect_oshape, top_face, rmesh3x2) == 0

@test Mesh.integral_face_rel_on_oshape_face(x^3*y^4, rect_oshape, Mesh.interior_face, rmesh3x2) == 1/20
@test Mesh.integral_face_rel_on_oshape_face(y^4, rect_oshape, left_face, rmesh3x2) == 1/5
@test Mesh.integral_face_rel_on_oshape_face(y^4, rect_oshape, right_face, rmesh3x2) == 1/5
@test Mesh.integral_face_rel_on_oshape_face(x*y^4, rect_oshape, right_face, rmesh3x2) == 0
@test Mesh.integral_face_rel_on_oshape_face(x^3, rect_oshape, bottom_face, rmesh3x2) == 1/4
@test Mesh.integral_face_rel_on_oshape_face(x^3, rect_oshape, top_face, rmesh3x2) == 1/4
@test Mesh.integral_face_rel_on_oshape_face(x^3*y, rect_oshape, top_face, rmesh3x2) == 0

## Integrate polynomial on faces of reference element.
@test nearly_eq(Mesh.integral_face_rel_on_oshape_face(2*x*y^2 + 3*x^3*y^4, rect_oshape, Mesh.interior_face, rmesh3x2), 2*(1/6) + 3*(1/20))
@test nearly_eq(Mesh.integral_face_rel_on_oshape_face(1.2*y^2 + 3.4*y^4, rect_oshape, left_face, rmesh3x2), 1.2/3 + 3.4/5)
@test nearly_eq(Mesh.integral_face_rel_on_oshape_face(1.5*y^2 + 0.123*y^4, rect_oshape, right_face, rmesh3x2), 1.5*(1/3) + 0.123*(1/5))
@test nearly_eq(Mesh.integral_face_rel_on_oshape_face(1.2x + 2x^3, rect_oshape, bottom_face, rmesh3x2), 1.2/2 + 2/4)
@test nearly_eq(Mesh.integral_face_rel_on_oshape_face(4.5*x*y^2 + 23.2x^3, rect_oshape, top_face, rmesh3x2), 23.2/4)

# Integrate vector monomials vs outward normals on the side faces of the reference finite element.
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(one_mon, VM(x^3*y^4, dim(1)), rect_oshape, left_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y, VM(x^3*y^4, dim(1)), rect_oshape, right_face, rmesh3x2), 1/6)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y, VM(y^2, dim(1)), rect_oshape, left_face, rmesh3x2), -1/4)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x*y, VM(y^2, dim(1)), rect_oshape, right_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y, VM(x*y^2, dim(1)), rect_oshape, bottom_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y, VM(y^2, dim(1)), rect_oshape, top_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x, VM(x^2, dim(2)), rect_oshape, bottom_face, rmesh3x2), -1/4)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x, VM(x^2, dim(2)), rect_oshape, top_face, rmesh3x2), 1/4)

# Test integration of a product of an arbitrary function and an element-local monomial on finite element faces.

f3(x::Vector{R}) = x[1]^2 * x[2]^3
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x*y, fenum(1), Mesh.interior_face, rmesh3x2), 1/20)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, y, fenum(1), left_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, y, fenum(1), right_face, rmesh3x2), 1/5)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x, fenum(1), bottom_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x, fenum(1), top_face, rmesh3x2), 1/4)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x*y, fenum(1), top_face, rmesh3x2), 0)

# Integrating the product below on fe 5 should be equivalent to integrating the monomial x^3 y^4 z on the reference element interior.
fe5_coords = Mesh.fe_interior_origin(fenum(5), rmesh3x2)
f4(x::Vector{R}) = (x[1] - fe5_coords[1])^2 * (x[2] - fe5_coords[2])^3
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x*y, fenum(5), Mesh.interior_face, rmesh3x2), 1/20)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, y, fenum(5), left_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, y, fenum(5), right_face, rmesh3x2), 1/5)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x, fenum(5), bottom_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x, fenum(5), top_face, rmesh3x2), 1/4)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x*y, fenum(5), top_face, rmesh3x2), 0)
*/

fn assert_approx(a:R, b:R) -> () {
  assert!(abs(a - b) < 10e-10)
}

fn main() {
  println("hello");
}

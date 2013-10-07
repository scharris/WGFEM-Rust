use std::vec;
use std::num::sqrt;
use std::iter::range_inclusive;

use common::*;
use monomial::{Monomial, Mon1d, Mon2d, Mon3d, Mon4d};
use polynomial::Polynomial;
use vector_monomial::VectorMonomial;
use mesh::*;
use quadrature::*;

// auxiliary types

// MeshCoord represents a single integer-valued mesh coordinates component,
// e.g. column or row or stack, etc.
#[deriving(Eq, TotalEq, Ord, TotalOrd, Clone)]
pub struct MeshCoord(uint);

pub struct NBSideGeom {
  perp_axis: Dim,
  mesh_coords: ~[MeshCoord]
}

pub struct RectMesh<Mon> {

  // The number of spatial dimensions of the Euclidiean space containing the mesh. 
  space_dims: uint,

  // Mesh coordinate ranges in R^d defining the boundaries of the mesh.
  min_bounds: ~[R],
  max_bounds: ~[R],

  // Logical dimensions of the mesh in integer mesh axis coordinates,
  // with directions corresponding to the coordinate axes (cols, rows,...).
  mesh_ldims: ~[MeshCoord],

  // Actual dimensions of any single finite element, the displacement vector from the
  // minimum coordinates corner to the maximum coordinates corner.
  fe_side_lens: ~[R],

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


fn new_impl<Mon:Monomial>(min_bounds: ~[R],
                          max_bounds: ~[R],
                          mesh_ldims: ~[MeshCoord],
                          integration_rel_err: R,
                          integration_abs_err: R) -> ~RectMesh<Mon> {

  let space_dims = Monomial::domain_space_dims(None::<Mon>); // TODO
  assert!(min_bounds.len() == space_dims);
  assert!(max_bounds.len() == space_dims);
  assert!(mesh_ldims.len() == space_dims);
  
  let fe_side_lens: ~[R] =
    vec::from_fn(space_dims, |r| {
      let bounds_diff = max_bounds[r] - min_bounds[r];
      let ldim_r = *mesh_ldims[r];
      assert!(bounds_diff > 0 as R);
      assert!(ldim_r > 0);
      bounds_diff/(ldim_r as R)
    });

  let cumprods_mesh_ldims: ~[uint] =
    mesh_ldims.iter().scan(1, |prod, &ldim| {
      *prod *= *ldim;
      Some(*prod)
    }).to_owned_vec();

  let cumprods_nb_side_mesh_ldims_by_perp_axis: ~[~[uint]] =
    vec::from_fn(space_dims, |perp_axis| {
      vec::from_fn(space_dims, |prods_top_dim| {
        range_inclusive(0, prods_top_dim).fold(1u, |ldims_prod, r| {
          ldims_prod * (if r != perp_axis { *mesh_ldims[r] } else { *mesh_ldims[r]-1 })
        })
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

  let rect_diameter = sqrt(fe_side_lens.iter().fold(0 as R, |sum_sq_lens, &len| sum_sq_lens + len*len));

  ~RectMesh {
    space_dims: space_dims,
    min_bounds: min_bounds,
    max_bounds: max_bounds,
    mesh_ldims: mesh_ldims,
    fe_side_lens: fe_side_lens,
    cumprods_mesh_ldims: cumprods_mesh_ldims,
    cumprods_nb_side_mesh_ldims_by_perp_axis: cumprods_nb_side_mesh_ldims_by_perp_axis,
    first_nb_side_nums_by_perp_axis: first_nb_side_nums_by_perp_axis,
    num_fes: num_fes,
    num_nb_sides: num_nb_sides,
    num_side_faces_per_fe: 2 * space_dims,
    rect_diameter: rect_diameter,
    rect_diameter_inv: 1./rect_diameter,
    integration_rel_err: integration_rel_err,
    integration_abs_err: integration_abs_err
  }
}


impl<Mon:Monomial> RectMesh<Mon> {
  
  pub fn new(min_bounds: ~[R],
             max_bounds: ~[R],
             mesh_ldims: ~[MeshCoord]) -> ~RectMesh<Mon> {
    new_impl(min_bounds, max_bounds, mesh_ldims,
             DEFAULT_INTEGRATION_REL_ERR, DEFAULT_INTEGRATION_ABS_ERR)
  }

  pub fn new_with_err_tols(min_bounds: ~[R],
                           max_bounds: ~[R],
                           mesh_ldims: ~[MeshCoord],
                           integration_rel_err: R,
                           integration_abs_err: R) -> ~RectMesh<Mon> {
      new_impl(min_bounds, max_bounds, mesh_ldims,
               integration_rel_err, integration_abs_err)
  }


  // side-related functions

  // Finds the perpendicular axis for a given non-boundary side in the mesh.
  #[inline]
  pub fn perp_axis_for_nb_side(&self, n: NBSideNum) -> Dim {
    assert!(*n < self.num_nb_sides);
    let mut r = self.space_dims-1;
    loop {
      if self.first_nb_side_nums_by_perp_axis[r] <= n { return Dim(r) }
      else if r == 0 { fail!("cannot find perpendicular axis for non-boundary side"); }
      else { r -= 1 };
    }
  }
  
  // Returns the geometric information for a non-boundary side in the mesh, which identifies the
  // perpendicular axis for the side, together with its coordinates in the orientation-specific mesh.
  pub fn nb_side_geom(&self, n: NBSideNum) -> NBSideGeom {
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
    let sidemesh_coords = vec::from_fn(self.space_dims, |r| {
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
  pub fn nb_side_with_mesh_coords(&self, coords: &[MeshCoord], perp_axis: Dim) -> NBSideNum {
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
    NBSideNum(range(1, self.space_dims).fold(*s_a0 + *coords[0], |sum_coord_contrs, r| {
      sum_coord_contrs + *coords[r] * self.cumprods_nb_side_mesh_ldims_by_perp_axis[*perp_axis][r-1]
    }))
  }

  // Converts finite element/interior coords in the main mesh to a finite element/interior number.
  #[inline]
  pub fn fe_with_mesh_coords(&self, coords: &[MeshCoord]) -> FENum {
    /*
     * The finite element (or interior) number for given mesh coordinates (c_1,...,c_d) is
     *   i_#(c_1,...,c_d) = sum_{i=1..d} { c_i prod_{l=1..i-1} k_l }
     *                   = c_1 + sum_{i=2..d} { (c_i - 1) prod_{l=1..i-1} k_l }
     * where k_l is the l^th component of the mesh dimensions.
     * Sum the contributions from the coordinates, each coordinate's value weighted with the
     * cumulative product of the logical mesh dimensions for lesser coordinate dimensions.
     * See Rectangular_Meshes.pdf document for the derivation.
     */
    FENum(range(1, self.space_dims).fold(*coords[0], |sum_coord_contrs, r| {
      sum_coord_contrs + *coords[r] * self.cumprods_mesh_ldims[r-1]
    }))
  }

  #[inline]
  pub fn fe_mesh_coords(&self, fe: FENum) -> ~[MeshCoord] {
    vec::from_fn(self.space_dims, |r| self.fe_mesh_coord(Dim(r), fe))
  }
  #[inline]
  pub fn fe_mesh_coord(&self, r: Dim, fe: FENum) -> MeshCoord {
    /*
     * The r^th 0-based mesh coordinate of side n is
     *   π(r,n) = (n mod (k_1 ··· k_r)) \ (k_1 ··· k_(r−1))
     * where k_i is the i^th component of the mesh dimensions.
     * See Rectangular_Meshes.pdf document for the derivation.
     */
    assert!(*r < self.space_dims);
    assert!(*fe < self.num_fes);
    let cumprods_preceeding_ldims = if *r == 0 { 1 } else { self.cumprods_mesh_ldims[*r-1] };
    MeshCoord((*fe % self.cumprods_mesh_ldims[*r]) / cumprods_preceeding_ldims)
  }
  
  #[inline]
  pub fn fe_coord_mins_corner(&self, fe: FENum) -> ~[R] {
    vec::from_fn(self.space_dims, |r| {
      let fe_mesh_coord_r = *self.fe_mesh_coord(Dim(r), fe) as R;
      self.min_bounds[r] + fe_mesh_coord_r * self.fe_side_lens[r]
    })
  }

}


impl<Mon:Monomial+RectIntegrable> Mesh<Mon>
                              for RectMesh<Mon> {

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
    assert!(*fe < self.num_fes);
    OShape(0)
  }
  
  #[inline(always)]
  fn num_side_faces_for_fe(&self, fe: FENum) -> uint {
    assert!(*fe < self.num_fes);
    self.num_side_faces_per_fe
  }
  
  #[inline(always)]
  fn num_side_faces_for_shape(&self, oshape: OShape) -> uint {
    assert!(*oshape == 0);
    self.num_side_faces_per_fe 
  }
  
  #[inline(always)]
  fn dependent_dim_for_oshape_side(&self, oshape: OShape, side_face: SideFace) -> Dim {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
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
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
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
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    let coord_a = self.fe_mesh_coord(a, fe);
    let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
    *coord_a == 0 && is_lesser_side || !is_lesser_side && *coord_a == *self.mesh_ldims[*a]-1
  }
  
  fn num_boundary_sides(&self) -> uint {
    range(0, self.space_dims).fold(0u, |perp_axis_contrs, perp_axis| {
      perp_axis_contrs +
      range(0, self.space_dims).fold(1u, |prod, r| {
        prod * if r == perp_axis { 2 } else { *self.mesh_ldims[r] }
      })
    })
  }
  
  #[inline(always)]
  fn shape_diameter_inv(&self, oshape: OShape) -> R {
    assert!(*oshape == 0);
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
    range(0, self.num_side_faces_for_fe(fe))
      .count(|sf| self.is_boundary_side(fe, SideFace(sf))) 
  }
  
  #[inline]
  fn max_num_shape_sides(&self) -> uint {
    self.num_side_faces_per_fe 
  }


  // integration functions

  #[inline]
  fn intg_global_fn_on_fe_int(&self, f: &fn(&[R]) -> R, fe: FENum) -> R {
    let d = self.space_dims;
    let fe_min_corner = self.fe_coord_mins_corner(fe);
    let fe_max_corner = vec::from_fn(d, |r| fe_min_corner[r] + self.fe_side_lens[r]);

    quadrature(&f, fe_min_corner, fe_max_corner, self.integration_rel_err, self.integration_abs_err)
  }

  #[inline]
  fn intg_global_fn_x_facerel_mon_on_fe_int(&self, f: &fn(&[R]) -> R, mon: Mon, fe: FENum) -> R {
    let d = self.space_dims;
    let fe_min_corner = &self.fe_coord_mins_corner(fe);
    let fe_int_origin = fe_min_corner;
    let fe_max_corner = vec::from_fn(d, |r| fe_min_corner[r] + self.fe_side_lens[r]);
    
    quadrature(&|x: &[R]| { f(x) * mon.value_at_for_origin(x, *fe_int_origin) },
               *fe_min_corner, fe_max_corner,
               self.integration_rel_err, self.integration_abs_err)
  }

  #[inline]
  fn intg_facerel_poly_on_oshape_int<P:Polynomial<Mon>>(&self, p: &P, oshape: OShape) -> R {
    assert!(*oshape == 0);
    p.foldl_terms(0 as R, |sum, (coef, mon)| {
      sum + coef * mon.integral_over_rect_at_origin(self.fe_side_lens)  
    })
  }

  #[inline]
  fn intg_facerel_poly_x_facerel_poly_on_oshape_int<P:Polynomial<Mon>>(&self, p1: &P, p2: &P, oshape: OShape) -> R {
    assert!(*oshape == 0);
    p1.foldl_terms(0 as R, |sum, (coef1, mon1)| {
      p2.foldl_terms(sum, |sum, (coef2, mon2)| {
        sum + coef1 * coef2 * (mon1*mon2).integral_over_rect_at_origin(self.fe_side_lens)
      })
    })
  }

  #[inline]
  fn intg_facerel_poly_x_facerel_poly_on_oshape_side<P:Polynomial<Mon>>(&self, p1: &P, p2: &P, oshape: OShape, side_face: SideFace) -> R {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    p1.foldl_terms(0 as R, |sum, (coef1, mon1)| {
      p2.foldl_terms(sum, |sum, (coef2, mon2)| {
        sum + coef1 * coef2 * (mon1*mon2).surface_integral_siderel_over_rect_side(self.fe_side_lens, a)
      })
    })
  }

  #[inline]
  fn intg_facerel_mon_on_oshape_int(&self, mon: Mon, oshape: OShape) -> R {
    assert!(*oshape == 0);
    mon.integral_over_rect_at_origin(self.fe_side_lens)
  }

  #[inline]
  fn intg_facerel_mon_on_oshape_side(&self, mon: Mon, oshape: OShape, side_face: SideFace) -> R {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    mon.surface_integral_siderel_over_rect_side(self.fe_side_lens, a)
  }

  #[inline]
  fn intg_facerel_mon_x_facerel_poly_on_oshape_int<P:Polynomial<Mon>>(&self, mon: Mon, p: &P, oshape: OShape) -> R {
    assert!(*oshape == 0);
    p.foldl_terms(0 as R, |sum, (coef, p_mon)| {
      sum + coef * (mon*p_mon).integral_over_rect_at_origin(self.fe_side_lens)
    })
  }

  #[inline]
  fn intg_facerel_mon_x_facerel_poly_on_oshape_side<P:Polynomial<Mon>>(&self, mon: Mon, p: &P, oshape: OShape, side_face: SideFace) -> R {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    p.foldl_terms(0 as R, |sum, (coef, p_mon)| {
      sum + coef * (mon*p_mon).surface_integral_siderel_over_rect_side(self.fe_side_lens, a)
    })
  }

  #[inline]
  fn intg_intrel_mon_x_siderel_mon_on_oshape_side(&self, int_mon: Mon, side_mon: Mon, oshape: OShape, side_face: SideFace) -> R {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
    let side_intrel_a_coord = if is_lesser_side { 0 as R } else { self.fe_side_lens[*a] };

    /* Here we break the interior-relative monomial on the side into the constant a-dim factor and the monomial
       of other dimension factors. Since the interior and side-relative coordinate systems differ only in dimension 
       a, the latter monomial has the same expression in side-relative coordinates. */
    let int_mon_dim_a_fac = pow(side_intrel_a_coord, *int_mon.exp(a) as uint);
    let int_mon_wo_dim_a_fac = int_mon.map_exp(a, |_| Deg(0));

    int_mon_dim_a_fac * (int_mon_wo_dim_a_fac * side_mon).surface_integral_siderel_over_rect_side(self.fe_side_lens, a)
  }

  #[inline]
  fn intg_siderel_mon_x_intrel_vmon_dot_normal_on_oshape_side(&self, side_mon: Mon, int_vmon: &VectorMonomial<Mon>, oshape: OShape, side_face: SideFace) -> R {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    match int_vmon.mon_dim() {
      Dim(r) if r == *a => {
        let int_vmon_mon = int_vmon.mon();
        let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
        let side_intrel_a_coord = if is_lesser_side { 0 as R } else { self.fe_side_lens[*a] };

        /* Here we break the interior-relative monomial on the side into the constant a-dim factor and the monomial
           of other dimension factors. Since the interior and side-relative coordinate systems differ only in dimension 
           a, the latter monomial has the same expression in side-relative coordinates. */
        let int_vmon_mon_dim_a_fac = pow(side_intrel_a_coord, *int_vmon_mon.exp(a) as uint);
        let int_vmon_mon_wo_dim_a_fac = int_vmon_mon.map_exp(a, |_| Deg(0));

        let outward_sense = if is_lesser_side { -1 as R } else { 1 as R };

        outward_sense * 
        int_vmon_mon_dim_a_fac *
        (int_vmon_mon_wo_dim_a_fac * side_mon).surface_integral_siderel_over_rect_side(self.fe_side_lens, a)
      }
      _ => 0 as R
    }
  }
 
  fn intg_siderel_poly_x_intrel_vmon_dot_normal_on_oshape_side<P:Polynomial<Mon>>(&self, p: &P, int_vmon: &VectorMonomial<Mon>, oshape: OShape, side_face: SideFace) -> R {
    assert!(*oshape == 0);
    assert!((*side_face as uint) < self.num_side_faces_per_fe);
    let a = side_face_perp_axis(side_face);
    match int_vmon.mon_dim() {
      Dim(r) if r == *a => {
        let int_vmon_mon = int_vmon.mon();
        let is_lesser_side = side_face_is_lesser_on_perp_axis(side_face);
        let side_intrel_a_coord = if is_lesser_side { 0 as R } else { self.fe_side_lens[*a] };

        /* Here we break the interior-relative monomial on the side into the constant a-dim factor and the monomial
           of other dimension factors. Since the interior and side-relative coordinate systems differ only in dimension 
           a, the latter monomial has the same expression in side-relative coordinates. */
        let int_vmon_mon_dim_a_fac = pow(side_intrel_a_coord, *int_vmon_mon.exp(a) as uint);
        let int_vmon_mon_wo_dim_a_fac = int_vmon_mon.map_exp(a, |_| Deg(0));

        let outward_sense = if is_lesser_side { -1 as R } else { 1 as R };

        outward_sense * 
        int_vmon_mon_dim_a_fac *
        p.foldl_terms(0 as R, |sum, (coef, mon)| {
          sum + coef * (int_vmon_mon_wo_dim_a_fac * mon).surface_integral_siderel_over_rect_side(self.fe_side_lens, a)
        })
      }
      _ => 0 as R
    }
  }

}


// RectIntegrable trait to be implemented by monomial types

pub trait RectIntegrable {

  /// Integrate the monomial over the rectangle of indicated dimensions having its minimums corner at the origin.
  fn integral_over_rect_at_origin(&self, side_lens: &[R]) -> R;  

  /// Integrate a *side-relative* monomial over a side face of a rectangle of given dimenions, the side being
  /// perpendicular to the indicated axis. The monomial's origin is the side's corner of minimum coordinates.
  /// Note that the integral will have the same value no matter which of the two sides perpendicular to the
  /// indicated axis is chosen as the domain of integration, because of the side local coordinates used to
  /// interpet the monomial (the  monomial will either have 0 exponent for the perpendicular axis variable
  /// and thus be independent of values in that coordinate, or else the integral will be 0).
  fn surface_integral_siderel_over_rect_side(&self, side_lens: &[R], side_perp_axis: Dim) -> R;

}

impl RectIntegrable for Mon1d {
  #[inline]
  fn integral_over_rect_at_origin(&self, side_lens: &[R]) -> R {
    let exp_plus_1 = *self.exps[0] as uint + 1;
    pow(side_lens[0], exp_plus_1)/(exp_plus_1 as R)
  }

  fn surface_integral_siderel_over_rect_side(&self, _: &[R], _: Dim) -> R {
    if *self.exps[0] != 0u8 {0 as R} else {1 as R}
  }
}

impl RectIntegrable for Mon2d {
  #[inline]
  fn integral_over_rect_at_origin(&self, side_lens: &[R]) -> R {
    let exp0_plus_1 = *self.exps[0] as uint + 1;
    let exp1_plus_1 = *self.exps[1] as uint + 1;
    pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
    pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R)
  }
  #[inline]
  fn surface_integral_siderel_over_rect_side(&self, side_lens: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[*side_perp_axis] != 0u8 { 0 as R }
    else {
      match side_perp_axis {
        Dim(0) => {
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R)
        }
        Dim(1) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R)
        }
        _ => fail!("Dimension out of range.")
      }
    }
  }
}

impl RectIntegrable for Mon3d {
  #[inline]
  fn integral_over_rect_at_origin(&self, side_lens: &[R]) -> R {
    let exp0_plus_1 = *self.exps[0] as uint + 1;
    let exp1_plus_1 = *self.exps[1] as uint + 1;
    let exp2_plus_1 = *self.exps[2] as uint + 1;

    pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
    pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R) *
    pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R)
  }
  #[inline]
  fn surface_integral_siderel_over_rect_side(&self, side_lens: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[*side_perp_axis] != 0u8 { 0 as R }
    else {
      match side_perp_axis {
        Dim(0) => {
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R)
        }
        Dim(1) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R)
        }
        Dim(2) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R)
        }
        _ => fail!("Dimension out of range.")
      }
    }
  }
}

impl RectIntegrable for Mon4d {
  #[inline]
  fn integral_over_rect_at_origin(&self, side_lens: &[R]) -> R {
    let exp0_plus_1 = *self.exps[0] as uint + 1;
    let exp1_plus_1 = *self.exps[1] as uint + 1;
    let exp2_plus_1 = *self.exps[2] as uint + 1;
    let exp3_plus_1 = *self.exps[3] as uint + 1;
    pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
    pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R) *
    pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R) *
    pow(side_lens[3], exp3_plus_1)/(exp3_plus_1 as R)
  }
  #[inline]
  fn surface_integral_siderel_over_rect_side(&self, side_lens: &[R], side_perp_axis: Dim) -> R {
    if *self.exps[*side_perp_axis] != 0u8 { 0 as R }
    else {
      match side_perp_axis {
        Dim(0) => {
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          let exp3_plus_1 = *self.exps[3] as uint + 1;
          pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R) *
          pow(side_lens[3], exp3_plus_1)/(exp3_plus_1 as R)
        }
        Dim(1) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          let exp3_plus_1 = *self.exps[3] as uint + 1;
          pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R) *
          pow(side_lens[3], exp3_plus_1)/(exp3_plus_1 as R)
        }
        Dim(2) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp3_plus_1 = *self.exps[3] as uint + 1;
          pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(side_lens[3], exp3_plus_1)/(exp3_plus_1 as R)
        }
        Dim(3) => {
          let exp0_plus_1 = *self.exps[0] as uint + 1;
          let exp1_plus_1 = *self.exps[1] as uint + 1;
          let exp2_plus_1 = *self.exps[2] as uint + 1;
          pow(side_lens[0], exp0_plus_1)/(exp0_plus_1 as R) *
          pow(side_lens[1], exp1_plus_1)/(exp1_plus_1 as R) *
          pow(side_lens[2], exp2_plus_1)/(exp2_plus_1 as R)
        }
        _ => fail!("Dimension out of range.")
      }
    }
  }
}

// side-related auxiliary stateless functions

// Find the axis which is perpendicular to the given side face.
#[inline]
fn side_face_perp_axis(side_face: SideFace) -> Dim {
  Dim((*side_face as uint) / 2)
}

// Determine whether a side face is the one with lesser axis value along its perpendicular axis.
#[inline(always)]
fn side_face_is_lesser_on_perp_axis(side_face: SideFace) -> bool {
  *side_face % 2 == 0
}

// Returns the side face of lesser coordinate value along the indicated axis.
#[inline]
pub fn lesser_side_face_perp_to_axis(a: Dim) -> SideFace {
  SideFace(2 * *a as uint)
}

// Returns the side face of greater coordinate value along the indicated axis.
#[inline]
pub fn greater_side_face_perp_to_axis(a: Dim) -> SideFace {
  SideFace((2 * *a + 1) as uint)
}


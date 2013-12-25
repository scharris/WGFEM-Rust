use common::*;
use monomial::{Monomial, Mon1d, Mon2d, Mon3d, Mon4d, domain_space_dims};
use polynomial::Polynomial;
use vector_monomial::VectorMonomial;
use mesh::*;
use quadrature::*;
use storage_by_ints::StorageByInts2;

use std::vec;
use std::num::{sqrt, hypot, max, abs, pow_with_uint};
use std::iter::{Iterator, range_inclusive, AdditiveIterator};
use std::cast;
use std::hashmap::{HashMap, HashSet};
use std::util::swap;


// Tag for a mesh geometric region or physical entity, from the input data.
#[deriving(Eq, TotalEq, Ord, TotalOrd, Clone)]
struct Tag(int);

// Node number from input data.
#[deriving(Eq, TotalEq, Ord, TotalOrd, Clone)]
struct NodeNum(uint);

type Point = (R,R);
type Vec = (R,R);



// Represents an element in the input which may be subdivided to contribute elements to the final mesh.
struct BaseTri {

  node_num_0: NodeNum,
  node_num_1: NodeNum,
  node_num_2: NodeNum,
  
  tag_phys_reg: Tag,
  tag_geom_ent: Tag,
  other_tags: Option<~[Tag]>

}

impl BaseTri {
  
  fn extra_subdiv_iters(&self) -> uint
  {
    // TODO: Allow specifying extra iterations via a Gmsh tag.
    0u
  }

  /* For each of the three sides of the indicated base mesh triangle to be subdivided, we can specify here the
     number of faces that the generated subtriangles should have between their vertexes which lie on this side.
     This allows these subdivision elements to meet those of a finer subdivision in a mesh element adjacent
     to this element ("hanging nodes").  The numbers are returned as a triplet of integers, corresponding to
     the face counts to be generated for elements' sides within sides v0v1, v1v2, and v2v0. */
  fn nums_side_faces_between_vertexes(&self) -> (u8,u8,u8)
  {
    // TODO: Allow specifying somewhere such as in mesh element tags.
    (1u8,1u8,1u8)
  }
  
}

// Reference triangle type representing the oriented shapes of the mesh.
struct RefTri {

  v01: Vec, // displacement vector from vertex 0 to vertex 1
  v02: Vec, // displacement vector from vertex 0 to vertex 2
 
  // There can be multiple side faces between vertexes, the numbers of which are recorded here for each vertex pair.
  nums_side_faces_between_vertexes: (u8,u8,u8), // v0 <-> v1, v1 <-> v2, v2 <-> v0
  
  num_side_faces: uint,
  outward_normals_by_side_face: ~[Vec],
  dep_dims_by_side_face: ~[Dim],

  diameter_inv: R,

}

// Constructor for reference triangles.
impl RefTri {

  pub fn new(v0: Point, v1: Point, v2: Point, scale: R, nums_side_faces_btw_verts: (u8,u8,u8)) -> RefTri
  {
    let ((v0_x,v0_y), (v1_x,v1_y), (v2_x,v2_y)) = (v0,v1,v2);
    let scaled_v01 = (scale*(v1_x - v0_x), scale*(v1_y - v0_y));
    let scaled_v02 = (scale*(v2_x - v0_x), scale*(v2_y - v0_y));
    let norm_scaled_v12 = scale * hypot(v2_x - v1_x, v2_y - v1_y);
    let normals = RefTri::side_face_outward_normals(v0, v1, v2, nums_side_faces_btw_verts);
    let dep_dims_by_side_face = RefTri::side_face_dep_dims(v0,v1,v2, nums_side_faces_btw_verts);
    let diameter_inv = 1./max(norm(scaled_v01), max(norm(scaled_v02), norm_scaled_v12));
    
    let num_side_faces = { let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts; (sfs_v01 + sfs_v12 + sfs_v20) as uint };
    
    RefTri{ v01: scaled_v01,
            v02: scaled_v02,
            nums_side_faces_between_vertexes: nums_side_faces_btw_verts,
            outward_normals_by_side_face: normals,
            dep_dims_by_side_face: dep_dims_by_side_face,
            diameter_inv: diameter_inv,
            num_side_faces: num_side_faces }
  }

  /*
   Outward Normals for Side Faces
   The outward normal for an inter-vertex vector w (extended to R^3 via 0 third component) is
     n = (w x (0,0,cc)) / |w|
       = cc (w_2, -w_1, 0) / |w|. // (using one based indexing on w)
   where cc = 1 if w is part of a clockwise traversal of the triangle, and -1 otherwise.
   For any pair of successive inter-vertex vectors v and w, cc can be computed as:
     cc = sgn(z(v x w)), where z is projection of component 3.
  */
  fn side_face_outward_normals(v0: Point, v1: Point, v2: Point, nums_side_faces_btw_verts: (u8,u8,u8)) -> ~[Vec]
  {
    let ((v0_x, v0_y), (v1_x, v1_y), (v2_x, v2_y)) = (v0, v1, v2);
    let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
    
    // inter-vertex vectors
    let v01 = (v1_x - v0_x, v1_y - v0_y);
    let v12 = (v2_x - v1_x, v2_y - v1_y);
    let v20 = (v0_x - v2_x, v0_y - v2_y);
    let ivvs = [(v01, sfs_v01),
                (v12, sfs_v12),
                (v20, sfs_v20)];
    
    let cc = if RefTri::oriented_counterclockwise(v01, v12) { 1. as R } else { -1. as R };
    
    let mut normals = vec::with_capacity((sfs_v01 + sfs_v12 + sfs_v20) as uint);
    for &((ivv_x, ivv_y), sfs) in ivvs.iter()
    {
      let len = hypot(ivv_x, ivv_y);
      let n = (cc * ivv_y/len, -cc * ivv_x/len);
      for sf in range(0, sfs as uint)
      {
        normals.push(n);
      }
    }
    normals
  }
  
  // Determine whether a counterclockwise rotation not exceeding 1/2 turn can transform the first vector to a positive multiple of the second.
  fn oriented_counterclockwise((u_x,u_y): Vec, (v_x,v_y): Vec) -> bool
  {
    u_x*v_y - u_y*v_x > 0. as R
  }


  fn side_face_dep_dims(v0: Point, v1: Point, v2: Point, nums_side_faces_btw_verts: (u8,u8,u8)) -> ~[Dim]
  {
    let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
    let vert_pairs = [((v0,v1), sfs_v01),
                      ((v1,v2), sfs_v12),
                      ((v2,v0), sfs_v20)];
    
    let mut dep_dims = vec::with_capacity((sfs_v01 + sfs_v12 + sfs_v20) as uint);

    for &((va,vb), sfs) in vert_pairs.iter()
    {
      let dep_dim = RefTri::side_face_dep_dim(va, vb);
      for i in range(0, sfs)
      {
        dep_dims.push(dep_dim);
      }
    }
    dep_dims
  }

  // Return a dependent dimension for the given side endpoints.
  fn side_face_dep_dim((va_x, va_y): Point, (vb_x, vb_y): Point) -> Dim
  {
    if abs(vb_x - va_x) >= abs(vb_y - va_y) { Dim(1) } else { Dim(0) }
  }

} // RefTri impl


// A single finite element in the final mesh.
struct ElTri {

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

  // number of oriented shapes
  num_oshapes: uint,

  // integration support members
  integrand_work_array: ~[R],
  space_dim_zeros: ~[R],
  space_dim_ones: ~[R],
  integration_rel_err: R,
  integration_abs_err: R,

  // element -> Tag maps
  phys_reg_tags_by_fenum: Option<~[Tag]>,
  geom_ent_tags_by_fenum: Option<~[Tag]>,

}



/****************************************************************************
 * Mesh construction.
 */

// Mesh Builder
struct MeshBuilder {
  
  oshapes: ~[RefTri],
  fes: ~[ElTri],

  side_reps_by_endpts: HashMap<(Point,Point),SideReps>,

  phys_reg_tags_by_fenum: Option<~[Tag]>,
  geom_ent_tags_by_fenum: Option<~[Tag]>,

  num_nb_sides: uint,
  num_b_sides: uint,

  // work buffers
  sfs_btw_verts_work_set: HashSet<(u8,u8,u8)>,
  sf_endpt_pairs_work_buf: Option<~[(Point,Point)]>,
}


impl MeshBuilder {
  
  pub fn make_mesh<Mon:Monomial,I:Iterator<BaseTri>>
         (base_pts_by_nodenum: ~[Point],
          mut base_tris_iter: &mut I,
          est_base_tris: uint,
          global_subdiv_iters: uint,
          integration_rel_err: R,
          integration_abs_err: R,
          load_phys_reg_tags: bool,
          load_geom_ent_tags: bool) -> TriMesh<Mon>
  {
    // Make estimates of the number of finite elements and reference triangles for storage allocation.
    // This estimate ignores optional additional subdivisions that may be specified for some input elements.
    let est_fes = est_base_tris * pow_with_uint(4, global_subdiv_iters);
    // We'll estimate 2 reference triangles for each mesh element if we are subdividing, 1 if not.
    // We ignore for this estimate any additional reference triangles for additional subdivisions.
    let est_ref_tris = (if global_subdiv_iters > 0 { 2 } else { 1 }) * est_base_tris;
    
    // The MeshBuilder structure holds data to be updated as base elements are processed.
    let mut bldr = MeshBuilder {
      fes: vec::with_capacity(est_fes),
      oshapes: vec::with_capacity(est_ref_tris),
      side_reps_by_endpts: HashMap::with_capacity((est_fes * 3)/2 as uint), // estimate assumes most fes have 3 sides, each in 2 fes
      phys_reg_tags_by_fenum: if load_phys_reg_tags { Some(vec::with_capacity(est_fes)) } else { None },
      geom_ent_tags_by_fenum: if load_geom_ent_tags { Some(vec::with_capacity(est_fes)) } else { None },
      num_nb_sides: 0,
      num_b_sides: 0,
      sfs_btw_verts_work_set: HashSet::with_capacity(20),
      sf_endpt_pairs_work_buf: Some(vec::with_capacity(6)),
    };
    
    // Process input mesh elements.
    for base_tri in base_tris_iter
    {
      bldr.process_base_tri(&base_tri, base_pts_by_nodenum, global_subdiv_iters);
    }

    //// Create the final non-boundary sides data structures based on the mapping of side endpoints to fe faces.
    let (nbsidenums_by_fe_face, nbsideincls_by_nbsidenum) = bldr.create_nb_sides_data();
    assert!(nbsideincls_by_nbsidenum.len() == bldr.num_nb_sides as uint);

    let num_fes = bldr.fes.len(); 
    let num_oshapes = bldr.oshapes.len();
    let space_dims = domain_space_dims::<Mon>();

    TriMesh {
      fes: bldr.take_fes(),
      oshapes: bldr.take_oshapes(),
      nbsidenums_by_fe_face: nbsidenums_by_fe_face,
      nbsideincls_by_nbsidenum: nbsideincls_by_nbsidenum,
      num_fes: num_fes,
      num_nb_sides: bldr.num_nb_sides as uint,
      num_b_sides: bldr.num_b_sides as uint,
      num_oshapes: num_oshapes,
      integrand_work_array: vec::from_elem(space_dims, 0 as R),
      space_dim_zeros: vec::from_elem(space_dims, 0 as R),
      space_dim_ones: vec::from_elem(space_dims, 1 as R),
      integration_rel_err: integration_rel_err,
      integration_abs_err: integration_abs_err,
      phys_reg_tags_by_fenum: bldr.take_phys_reg_tags(),
      geom_ent_tags_by_fenum: bldr.take_geom_ent_tags(),
    }
  }

  fn process_base_tri(& mut self,
                      base_tri: &BaseTri,
                      mesh_pts_by_nodenum: &[Point],
                      global_subdiv_iters: uint)
  {
    let (v0, v1, v2) = (mesh_pts_by_nodenum[*base_tri.node_num_0],
                        mesh_pts_by_nodenum[*base_tri.node_num_1],
                        mesh_pts_by_nodenum[*base_tri.node_num_2]);

    // For each of the three sides of this mesh element to be subdivided, the generated subdivision
    // elements can be made to have more than one face between their vertexes lying on the side.  This is
    // useful to support "hanging" nodes where a finer subdivision is adjacent to this one.  The triplet
    // returned represents the number of side faces that a generated element has if its vertexes lie between
    // base triangle vertex pairs v0 and v1, v1 and v2, and v2 and v0, respectively.
    let nums_sfs_btw_verts = base_tri.nums_side_faces_between_vertexes();
    
    // Add any extra subdivision iterations to be done in this element.
    let subdiv_iters = global_subdiv_iters + base_tri.extra_subdiv_iters();

    // Register the primary reference triangles for our mesh element's subdivisions. Normally there will only be one
    // such primary (ie. non-inverted) reference triangle, however when multiple side faces are present between two vertexes
    // then additional reference triangles are required.
    let first_new_pri_oshape = self.register_primary_ref_tris(v0,v1,v2, nums_sfs_btw_verts, subdiv_iters);
    let last_new_pri_oshape = OShape(self.oshapes.len()-1);

    // Create a lookup function for the new primary oshape numbers by their numbers of side faces between vertexes.
    let pri_oshapes_by_nums_sfs_btw_verts: |(u8,u8,u8)| -> OShape = |sfs_btw_verts: (u8,u8,u8)| {
      for os in range_inclusive(*first_new_pri_oshape, *last_new_pri_oshape)
      {
        if self.oshapes[os].nums_side_faces_between_vertexes == sfs_btw_verts { return OShape(os); }
      }
      fail!("Reference triangle not found by numbers of side faces between vertexes.");
    };

    if subdiv_iters > 0
    {
      // We're subdividing, so register the secondary reference triangle.
      let sec_oshape = self.register_secondary_ref_tri(v0,v1,v2, subdiv_iters);

      // Do the subdivisions.
      self.subdivide_primary(v0, v1, v2,
                             subdiv_iters,
                             nums_sfs_btw_verts,
                             pri_oshapes_by_nums_sfs_btw_verts, sec_oshape,
                             base_tri.tag_phys_reg, base_tri.tag_geom_ent);
    }
    else // no subdivision to be done
    { 
      // The base triangle itself is our finite element, with its own reference triangle (added in register_primary_ref_tris above).
      self.add_fe(first_new_pri_oshape, v0,v1,v2, base_tri.tag_phys_reg, base_tri.tag_geom_ent);
    }
  }

  /* Create primary reference triangles for the given triangle to be subdivided. These reference triangles
     are rescaled translations of the original undivided triangle. If nums_side_faces_btw_verts is other than
     (1,1,1), then multiple reference triangles will be generated, differing in the number of side faces
     between their vertexes for support of "hanging" nodes.  The function returns a function mapping the
     numbers of side faces between vertexes to the newly registered reference triangle (oshape) number.
     Pains are taken to not create more reference triangles than are actually used, so that some global
     mesh properties such as maximum element diameter can be determined from only the reference elements.
     Returns the first new oriented shape number that was created. */
  fn register_primary_ref_tris(& mut self, v0: Point, v1: Point, v2: Point, nums_sfs_btw_verts: (u8,u8,u8), subdiv_iters: uint) -> OShape {

    let first_new_oshape = OShape(self.oshapes.len());
   
    let scale = { let p: uint = pow_with_uint(2u, subdiv_iters); 1./(p as R) };

    if subdiv_iters == 0 || nums_sfs_btw_verts == (1u8, 1u8, 1u8) 
    {
      self.oshapes.push(RefTri::new(v0,v1,v2, scale, nums_sfs_btw_verts));
    } 
    else // We are subdividing, and there is more than one side face between some vertex pair.
    {
      let (nums_sfs_btw_verts_0, nums_sfs_btw_verts_1, nums_sfs_btw_verts_2) = nums_sfs_btw_verts;

      // We need to add a separate primary reference triangle for each triplet representing the number of side
      // faces between vertexes that occur in primary subdivision triangles of this base triangle.
      // We use the set below to track the unique triplets of numbers of side faces between vertexes.
      self.sfs_btw_verts_work_set.clear();
    
      // Handle subdivision triangles at one of the base triangle's corners, v0, v1 or v2.
      // To support the corner elements, we need primary reference triangles which get 2 of their 3 numbers
      // of faces between vertex pairs from nums_side_faces_btw_verts, because they have two sides on the
      // original undivided triangle's sides.
      self.sfs_btw_verts_work_set.insert((nums_sfs_btw_verts_0, 1u8, nums_sfs_btw_verts_2)); // v0 corner element
      self.sfs_btw_verts_work_set.insert((nums_sfs_btw_verts_0, nums_sfs_btw_verts_1, 1u8)); // v1 corner element
      self.sfs_btw_verts_work_set.insert((1u8, nums_sfs_btw_verts_1, nums_sfs_btw_verts_2)); // v2 corner element

      // If subdividing more than once, then we also need:
      //  - a primary reference triangle with only one side face between each of its vertex pairs (within the central secondary triangle).
      //  - primary reference triangles which only have one pair of vertexes contained in one of the original base triangle's sides,
      //    so that one pair of vertexes may have multiple side faces between them but not so for the other two vertex pairs.
      if subdiv_iters > 1
      {
        self.sfs_btw_verts_work_set.insert((1u8, 1u8, 1u8));
        self.sfs_btw_verts_work_set.insert((nums_sfs_btw_verts_0, 1u8, 1u8));
        self.sfs_btw_verts_work_set.insert((1u8, nums_sfs_btw_verts_1, 1u8));
        self.sfs_btw_verts_work_set.insert((1u8, 1u8, nums_sfs_btw_verts_2));
      } 

      // Now create the reference triangles for each unique triplet of numbers of side faces between vertexes.
      for &sfs_btw_verts in self.sfs_btw_verts_work_set.iter()
      {
        self.oshapes.push(RefTri::new(v0,v1,v2, scale, sfs_btw_verts));
      }
    }
    
    first_new_oshape 
  }

  // Register the single secondary (inverted) subdivision triangle for the given (primary) base triangle, and number of subdivisions
  // to be applied on the base (primary) triangle.
  fn register_secondary_ref_tri(& mut self,
                                (pv0_x, pv0_y): Point, (pv1_x, pv1_y): Point, (pv2_x, pv2_y): Point, // base (primary) triangle vertexes
                                subdiv_iters: uint) -> OShape
  {
    let sec_oshape = OShape(self.oshapes.len());
    let sec_ref_tri = {
      let (sv0, sv1, sv2) = ((0.5*(pv0_x + pv1_x), 0.5*(pv0_y + pv1_y)),
                             (0.5*(pv1_x + pv2_x), 0.5*(pv1_y + pv2_y)),
                             (0.5*(pv2_x + pv0_x), 0.5*(pv2_y + pv0_y)));
      let scale = { let p: uint = pow_with_uint(2u, subdiv_iters-1); 1./(p as R) }; // sv0, sv1, sv2 already represents one subdivision
      RefTri::new(sv0,sv1,sv2, scale, (1u8,1u8,1u8))
    };
    self.oshapes.push(sec_ref_tri);
    sec_oshape
  }

  fn subdivide_primary(& mut self,
                       v0: Point, v1: Point, v2: Point,
                       iters: uint,
                       nums_side_faces_btw_verts: (u8,u8,u8),
                       pri_oshapes_by_nums_sfs_btw_verts: |(u8,u8,u8)| -> OShape,
                       sec_oshape: OShape,
                       tag_phys_reg: Tag, tag_geom_ent: Tag)
  {
    if iters == 0
    {
      let oshape = pri_oshapes_by_nums_sfs_btw_verts(nums_side_faces_btw_verts);
      self.add_fe(oshape, v0,v1,v2, tag_phys_reg, tag_geom_ent);
    }
    else
    {
      let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;

      let ((v0_x,v0_y), (v1_x,v1_y), (v2_x,v2_y)) = (v0,v1,v2);
      let midpt_v01 = (0.5*(v0_x + v1_x), 0.5*(v0_y + v1_y)); 
      let midpt_v12 = (0.5*(v1_x + v2_x), 0.5*(v1_y + v2_y));
      let midpt_v20 = (0.5*(v2_x + v0_x), 0.5*(v2_y + v0_y)); 


      // The sub-triangle including v0 has its first and last vertex pairs embedded in the original triangle's first and
      // third sides, and so inherits the corresponding numbers of faces between vertexes from nums_side_faces_btw_verts.
      // The middle vertex pair (midpt_v01 to midpt_v20) of this sub-triangle does not lie along an original side, however,
      // so has only one face between these vertexes.  Similar arguments apply for the remaining sub-triangles.
      self.subdivide_primary(v0, midpt_v01, midpt_v20,
                             iters-1,
                             (sfs_v01, 1, sfs_v20),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag_phys_reg, tag_geom_ent);

      self.subdivide_primary(midpt_v01, v1, midpt_v12,
                             iters-1,
                             (sfs_v01, sfs_v12, 1),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag_phys_reg, tag_geom_ent);

      self.subdivide_primary(midpt_v20, midpt_v12, v2,
                             iters-1,
                             (1, sfs_v12, sfs_v20),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag_phys_reg, tag_geom_ent);

      // secondary sub-triangle
      self.subdivide_secondary(midpt_v01, midpt_v12, midpt_v20,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag_phys_reg, tag_geom_ent);
    }
  }

  fn subdivide_secondary(& mut self,
                         v0: Point, v1: Point, v2: Point, // secondary triangle vertexes
                         iters: uint,
                         pri_oshapes_by_nums_sfs_btw_verts: |(u8,u8,u8)| -> OShape,
                         sec_oshape: OShape,
                         tag_phys_reg: Tag, tag_geom_ent: Tag)
  {
    if iters == 0
    { 
      self.add_fe(sec_oshape, v0,v1,v2, tag_phys_reg, tag_geom_ent);
    }
    else
    {
      let ((v0_x,v0_y), (v1_x,v1_y), (v2_x,v2_y)) = (v0,v1,v2);
      let midpt_v01 = (0.5*(v0_x + v1_x), 0.5*(v0_y + v1_y)); 
      let midpt_v12 = (0.5*(v1_x + v2_x), 0.5*(v1_y + v2_y));
      let midpt_v20 = (0.5*(v2_x + v0_x), 0.5*(v2_y + v0_y)); 

      // Corner subtriangles of this secondary triangle are secondary triangles needing one less subdivision iteration.
      self.subdivide_secondary(v0, midpt_v01, midpt_v20,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag_phys_reg, tag_geom_ent);

      self.subdivide_secondary(midpt_v01, v1, midpt_v12,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag_phys_reg, tag_geom_ent);

      self.subdivide_secondary(midpt_v20, midpt_v12, v2,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag_phys_reg, tag_geom_ent);

      // The central sub-triangle is a primary triangle.  We need to order the vertexes properly, so that all primary and
      // secondary triangles will have the same vertex orientations, in the sense that for any two subtriangles of the same
      // type (primary or secondary), some translation can take the vertexes of one to the corresponding vertexes of the other.
      // We do this by numbering the primary vertexes just as they were numbered in the original undivided primary triangle,
      // and numbering all secondary vertexes just as in the original central secondary triangle of the undivided primary
      // triangle.
      // We shifted "forward" when numbering vertexes of a central secondary subtriangle from the enclosing primary triangle,
      // so that e.g. the secondary's v0 is the midpoint of enclosing primary's v0 and v1. So now we must shift "back" for
      // this enclosed primary triangle, so e.g. the central primary triangle's v0 is the midpoint of this enclosing secondary
      // triangle's v2 and v0, and similarly for the other vertexes.
      self.subdivide_primary(midpt_v20, midpt_v01, midpt_v12,
                             iters-1,
                             (1u8,1u8,1u8),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag_phys_reg, tag_geom_ent);
    }
  }

  // This procedure will be used to create all finite elements during the processing of the input base triangles.
  fn add_fe(& mut self, oshape: OShape, v0: Point, v1: Point, v2: Point, tag_phys_reg: Tag, tag_geom_ent: Tag)
  {
    // Create the finite element.
    let fe = FENum(self.fes.len());
    self.fes.push(ElTri{oshape: oshape, v0: v0});
    
    // optionally store tags
    match self.phys_reg_tags_by_fenum { Some(ref mut tags) => tags.push(tag_phys_reg), _ => {} }
    match self.geom_ent_tags_by_fenum { Some(ref mut tags) => tags.push(tag_geom_ent), _ => {} }

    // Register this fe's local side representations ( fe/sf pairs ) by their endpoint pairs.
    let (nb_sides_delta, b_sides_delta) = {
      let nums_side_faces_btw_verts = self.oshapes[*oshape].nums_side_faces_between_vertexes;
      self.register_fe_side_reps_by_endpoints(fe, v0,v1,v2, nums_side_faces_btw_verts)
    };

    self.num_nb_sides = ((self.num_nb_sides as int) + nb_sides_delta) as uint;
    self.num_b_sides =  ((self.num_b_sides as int)  + b_sides_delta) as uint;
  }
  
  // Register the given finite element's side faces by their endpoints in the passed registry.
  fn register_fe_side_reps_by_endpoints(& mut self, 
                                        fe: FENum,
                                        v0: Point, v1: Point, v2: Point,
                                        nums_side_faces_btw_verts: (u8,u8,u8)) -> (int, int)
  {
    let mut nb_sides_delta = 0;
    let mut b_sides_delta = 0;
    
    // Take the endpoint pairs work buffer, filled with the side face endpoint pairs.
    let endpt_pairs = MeshBuilder::fill_side_face_endpoint_pairs(self.take_side_face_endpoint_pairs_buf(),
                                                                 v0,v1,v2, nums_side_faces_btw_verts, LesserEndpointsFirst);

    // Register the side face representations for this element by their endpoints.
    for sf in range(0, endpt_pairs.len())
    {
      let new_rep = (fe, SideFace(sf));
      self.side_reps_by_endpts.mangle(endpt_pairs[sf], (/*no context value*/),
        |_k, _ctx| { // no existing side rep found: insert
          b_sides_delta += 1;
          SideReps::new(new_rep) // new value to insert
        },
        |_k, side_reps, _ctx| { // existing side rep found: mutate
          assert!(side_reps.snd_rep.is_none());
          side_reps.add(new_rep);
          b_sides_delta -= 1;
          nb_sides_delta += 1;
        }
      );
    }

    // Give back ownership of the endpoint pairs work buffer.
    self.return_side_face_endpoint_pairs_buf(endpt_pairs);
    
    (nb_sides_delta, b_sides_delta)
  }

  /* This function defines the side faces enumeration for a finite element of given vertexes and numbers
     of faces between vertexes. Side faces are returned as an array of side endpoint pairs indexed by
     side face number. If lesser_endpts_first is true, then each endpoint pair endpoint will have
     the lesser point (compared lexicographically) in the first component of the pair. */
  fn fill_side_face_endpoint_pairs(mut endpts_buf: ~[(Point,Point)],
                                   v0: Point, v1: Point, v2: Point,
                                   nums_side_faces_btw_verts: (u8,u8,u8),
                                   endpoints_ordering: SideEndpointsOrdering) -> ~[(Point,Point)]
  {
    endpts_buf.clear();

    let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
    let side_faces = [(v0, v1, sfs_v01),
                      (v1, v2, sfs_v12),
                      (v2, v0, sfs_v20)];

    for &(va, vb, num_sfs_btw_va_vb) in side_faces.iter()
    {
      match num_sfs_btw_va_vb
      {
        1 => endpts_buf.push(mk_endpoint_pair(va, vb, endpoints_ordering)),
        2 => {
          let midpt = { let ((va_x,va_y),(vb_x,vb_y)) = (va,vb); (0.5*(va_x + vb_x), 0.5*(va_y + vb_y)) };
          endpts_buf.push(mk_endpoint_pair(va, midpt, endpoints_ordering));
          endpts_buf.push(mk_endpoint_pair(midpt, vb, endpoints_ordering));
        }
        _ => fail!("Only 1 or 2 faces between triangle vertexes are currently supported.")
      }
    }

    endpts_buf
  }

  // Based on the fe/side-face representations by endpoints (ignoring the endpoints), return:
  //   1) nb side numbers by (fe,face), newly assigned here
  //   2) an array of NBSideInclusions indexed by nb side number.
  fn create_nb_sides_data(& mut self) -> (StorageByInts2<Option<NBSideNum>>, ~[NBSideInclusions])
  {
    let max_sides_per_fe = self.oshapes.iter().map(|rt| rt.num_side_faces).max().unwrap();
    
    let mut nbsidenums_by_fe_face = StorageByInts2::from_elem(self.fes.len(), max_sides_per_fe, None);
    let mut nbsideincls_by_nbsidenum = vec::with_capacity(self.num_nb_sides);

    for (_, side_reps) in self.side_reps_by_endpts.iter()
    {
      match side_reps.snd_rep
      {
        Some((fe_2, sf_2)) => { // two reps: non-boundary side
          let (fe_1, sf_1) = side_reps.fst_rep;
          
          // Assign a new non-boundary side number.
          let nb_side_num = NBSideNum(nbsideincls_by_nbsidenum.len());

          // Add nb side inclusions structure for this nb side number.
          nbsideincls_by_nbsidenum.push(NBSideInclusions::new(nb_side_num, fe_1, sf_1, fe_2, sf_2));

          // Map the two fe/face pairs to this nb side number.
          nbsidenums_by_fe_face.set(*fe_1, *sf_1, Some(nb_side_num));
          nbsidenums_by_fe_face.set(*fe_2, *sf_2, Some(nb_side_num));
        }
        None => {} // boundary side: we're only interested in non-boundary sides
      }
    }

    (nbsidenums_by_fe_face, nbsideincls_by_nbsidenum)
  }


  // Functions for taking completed data from the builder, and borrowing/returning work buffers.

  fn take_fes(& mut self) -> ~[ElTri]
  {
    let mut fes = ~[];
    swap(& mut fes, & mut self.fes);
    fes
  }
  
  fn take_oshapes(& mut self) -> ~[RefTri]
  {
    let mut oshapes = ~[];
    swap(& mut oshapes, & mut self.oshapes);
    oshapes
  }
  
  fn take_phys_reg_tags(& mut self) -> Option<~[Tag]>
  {
    let mut tags = None;
    swap(& mut tags, & mut self.phys_reg_tags_by_fenum);
    tags 
  }
  
  fn take_geom_ent_tags(& mut self) -> Option<~[Tag]>
  {
    let mut tags = None;
    swap(& mut tags, & mut self.geom_ent_tags_by_fenum);
    tags 
  }

  fn take_side_face_endpoint_pairs_buf(& mut self) -> ~[(Point,Point)]
  {
    let mut endpt_pairs_opt = None;
    swap(& mut endpt_pairs_opt, & mut self.sf_endpt_pairs_work_buf);
    let mut endpt_pairs = endpt_pairs_opt.unwrap();
    endpt_pairs
  }
  
  fn return_side_face_endpoint_pairs_buf(& mut self, buf: ~[(Point,Point)])
  {
    let mut ret_buf_opt = Some(buf);
    swap(& mut ret_buf_opt, & mut self.sf_endpt_pairs_work_buf);
  }

} // MeshBuilder impl


// A structure to hold the one or two local representations of any side as fe/side face pairs, for use during mesh construction.
struct SideReps {
  fst_rep: (FENum, SideFace),
  snd_rep: Option<(FENum, SideFace)>,
}

impl SideReps {
  
  fn new(fe_sf: (FENum,SideFace)) -> SideReps { SideReps { fst_rep: fe_sf, snd_rep: None } } 

  fn add(& mut self, (fe,sf): (FENum,SideFace))
  {
    // Keep the representations in fe order.
    match self.fst_rep
    {
      (fe1,_) if fe1 < fe => { self.snd_rep = Some((fe,sf)); }
      (fe1,_) if fe1 > fe => { 
        self.snd_rep = Some(self.fst_rep);
        self.fst_rep = (fe,sf);
      }
      _ => fail!("Encountered side with two representations in the same finite element.")
    }
  }

}

/*
# construction from rectangle mesh
function from_rect_mesh(rmesh::RMesh.RectMesh, subdiv_ops::Int)
  assert(Mesh.space_dim(rmesh) == 2)
  const EMPTY_TAGS_ARRAY = Array(Tag,0)

  const num_rects = Mesh.num_fes(rmesh)
  const rect_dims = RMesh.fe_dims(rmesh)
  const num_points = uint64(num_rects + sum(RMesh.logical_dims(rmesh))+1)

  const nodenums_by_mcoords = sizehint(Dict{(RMesh.MeshCoord, RMesh.MeshCoord), NodeNum}(), num_points)
  const points_by_ptnum = sizehint(Array(Point, 0), num_points)

  function register_point(pt_mcoords::(RMesh.MeshCoord, RMesh.MeshCoord))
    const existing_ptnum = get(nodenums_by_mcoords, pt_mcoords, nodenum(0))
    if existing_ptnum == 0
      # Register the new point number with its mesh coordinates.
      push!(points_by_ptnum,
            Point(rmesh.min_bounds[1] + (pt_mcoords[1] - 1) * rect_dims[1],
                  rmesh.min_bounds[2] + (pt_mcoords[2] - 1) * rect_dims[2]))
      const new_ptnum = nodenum(length(points_by_ptnum))
      nodenums_by_mcoords[pt_mcoords] = new_ptnum
      new_ptnum
    else
      existing_ptnum
    end
  end

  const base_tris = sizehint(Array(BaseTri, 0), 2*num_rects)
  
  const rect_mcoords = Array(RMesh.MeshCoord, Mesh.space_dim(rmesh))

  for rect_fe=Mesh.fenum(1):num_rects
    RMesh.fe_mesh_coords!(rect_fe, rect_mcoords, rmesh)

    const ll_ptnum = register_point((rect_mcoords[1],   rect_mcoords[2]))
    const lr_ptnum = register_point((rect_mcoords[1]+1, rect_mcoords[2]))
    const ur_ptnum = register_point((rect_mcoords[1]+1, rect_mcoords[2]+1))
    const ul_ptnum = register_point((rect_mcoords[1],   rect_mcoords[2]+1))

    push!(base_tris, BaseTri((ll_ptnum, lr_ptnum, ur_ptnum), tag(rect_fe), tag(rect_fe), EMPTY_TAGS_ARRAY))
    push!(base_tris, BaseTri((ll_ptnum, ur_ptnum, ul_ptnum), tag(rect_fe), tag(rect_fe), EMPTY_TAGS_ARRAY))
  end

  TriMesh(points_by_ptnum,
          base_tris,
          length(base_tris),
          subdiv_ops,
          rmesh.integration_rel_err, rmesh.integration_abs_err,
          false, true) # load geom entity tags
end


# Mesh Construction
####################################################################



####################################################################
# implementations of required AbstractMesh functions

import Mesh.space_dim
space_dim(mesh::TriMesh) = SPACE_DIM

import Mesh.one_mon
one_mon(mesh::TriMesh) = mesh.one_mon

import Mesh.num_fes
num_fes(mesh::TriMesh) = mesh.num_fes

import Mesh.num_nb_sides
num_nb_sides(mesh::TriMesh) = mesh.num_nb_sides

import Mesh.num_oriented_element_shapes
num_oriented_element_shapes(mesh::TriMesh) = mesh.num_oshapes

import Mesh.oriented_shape_for_fe
oriented_shape_for_fe(fe::FENum, mesh::TriMesh) = mesh.fes[fe].oshape

import Mesh.num_side_faces_for_fe
num_side_faces_for_fe(fe::FENum, mesh::TriMesh) = ref_tri_for_fe(fe, mesh).num_side_faces

import Mesh.num_side_faces_for_shape
num_side_faces_for_shape(oshape::OShape, mesh::TriMesh) = mesh.oshapes[oshape].num_side_faces

import Mesh.dependent_dim_for_oshape_side
dependent_dim_for_oshape_side(oshape::OShape, sf::FEFaceNum, mesh::TriMesh) =
  mesh.oshapes[oshape].dep_dims_by_side_face[sf]

import Mesh.fe_inclusions_of_nb_side
fe_inclusions_of_nb_side(nbsn::NBSideNum, mesh::TriMesh) = mesh.nbsideincls_by_nbsidenum[nbsn]

import Mesh.nb_side_num_for_fe_side
nb_side_num_for_fe_side(fe::FENum, sf::FEFaceNum, mesh::TriMesh) = mesh.nbsidenums_by_fe_face[fe, sf]

import Mesh.is_boundary_side
is_boundary_side(fe::FENum, face::FEFaceNum, mesh::TriMesh) = mesh.nbsidenums_by_fe_face[fe, face] == 0

import Mesh.num_boundary_sides
num_boundary_sides(mesh::TriMesh) = mesh.num_b_sides

import Mesh.shape_diameter_inv
shape_diameter_inv(oshape::OShape, mesh::TriMesh) = mesh.oshapes[oshape].diameter_inv

import Mesh.max_fe_diameter
max_fe_diameter(mesh::TriMesh) = mapreduce(rt -> 1./rt.diameter_inv, max, mesh.oshapes)

import Mesh.fe_interior_origin!
function fe_interior_origin!(fe::FENum, fill::Vector{R}, mesh::TriMesh)
  const el_tri = mesh.fes[fe]
  fill[1] = el_tri.v1._1
  fill[2] = el_tri.v1._2
end


# Integration Functions

# Integrate a monomial on the indicated face of a finite element of given oriented
# shape, with the monomial interpreted locally on the face.
import Mesh.integral_face_rel_on_oshape_face
function integral_face_rel_on_oshape_face(mon::Monomial,
                                          oshape::OShape, face::FEFaceNum,
                                          mesh::TriMesh)
  const ref_tri = mesh.oshapes[oshape]

  if face == Mesh.interior_face
    # Order the interior-relative vertex points by their first coordinate values.
    const pts = sort!([zeroPt, ref_tri.v12, ref_tri.v13])
    const slope_12, slope_13 = slope_between(pts[1], pts[2]), slope_between(pts[1], pts[3])

    if pts[1]._1 == pts[2]._1 # vertical side on left between points 1 and 2
      # Integrate over the area bounded above and below by the lines between points 1 and 3 and 2 and 3,
      # and horizontally between the vertical left side formed by points 1 and 2, and point 3 on the right
      # where the upper and lower bounding lines meet.
      const slope_23 = slope_between(pts[2], pts[3])
      integral_mon_between_lines_meeting_at_point_and_vert_line(mon, pts[3], slope_13, slope_23, pts[1]._1)
    else
      # Points 1 and 2 do not form a vertical line. Integrate between points 1 and 2, and between points 2 and 3 if
      # points 2 and 3 don't lie on a vertical line.
      const fst_seg = integral_mon_between_lines_meeting_at_point_and_vert_line(mon, pts[1], slope_12, slope_13, pts[2]._1)
      const snd_seg =
        if pts[2]._1 == pts[3]._1 # vertical right side, no second segment
          zeroR
        else
          const slope_23 = slope_between(pts[3], pts[2])
          integral_mon_between_lines_meeting_at_point_and_vert_line(mon, pts[3], slope_13, slope_23, pts[2]._1)
        end
      fst_seg + snd_seg
    end
  else # line integral along side face
    # We want to compute int_{t=0..1} mon(p(t)-o) |p'(t)| dt, where p:[0,1] -> R^2
    # is a bijection traversing the side face smoothly, p(0) and p(1) being the side
    # endpoints, and o is the local origin for the side, which is the side's midpoint.
    const a,b = side_face_endpoint_pair(face,
                                       zeroPt, ref_tri.v12, ref_tri.v13,
                                       ref_tri.nums_side_faces_between_vertexes,
                                       false) # lesser endpoints first => false
    # The local origin of each side face is the midpoint.
    const o = 1/2*(a + b)
    # Path traversing the side from a to b, in local coordinates, relative to the side's local origin.
    const path_orel = let t = MON_VAR_1D
      [a._1 + (b._1 - a._1)*t - o._1,
       a._2 + (b._2 - a._2)*t - o._2]
    end
    const pullback_poly_1d = Poly.precompose_with_poly_path(mon, path_orel) * norm(b-a)
    const pullback_antider = Poly.antideriv(dim(1), pullback_poly_1d)
    Poly.value_at(pullback_antider, oneR) - Poly.value_at(pullback_antider, zeroR)
  end
end


# Integrate a global function f on the indicated finite element face, multiplied against
# a monomial interpreted locally on the face.
import Mesh.integral_global_x_face_rel_on_fe_face
function integral_global_x_face_rel_on_fe_face(f::Function,
                                               mon::Monomial,
                                               fe::FENum, face::FEFaceNum,
                                               mesh::TriMesh)
  const ref_tri = ref_tri_for_fe(fe, mesh)
  const fe_v1 = mesh.fes[fe].v1

  if face == Mesh.interior_face

    function f_x_mon(int_rel_x::Vector{R})
      const irel_1, irel_2 = int_rel_x[1], int_rel_x[2]
      const mon_val = Poly.value_at(mon, int_rel_x)
      try
        const x = int_rel_x # temporarily reuse storage of int_rel_x for global x
        # vertex 1 is the local origin for the interior
        x[1] += fe_v1._1
        x[2] += fe_v1._2
        f(x) * mon_val
      finally
        # restore int_rel_x to original contents
        int_rel_x[1] = irel_1; int_rel_x[2] = irel_2
      end
    end

    # Order the interior-relative vertex points by their first coordinate values.
    const pts = sort!([zeroPt, ref_tri.v12, ref_tri.v13])
    const slope_12, slope_13 = slope_between(pts[1], pts[2]), slope_between(pts[1], pts[3])

    if pts[1]._1 == pts[2]._1 # vertical side on left between points 1 and 2
      # Integrate over the area bounded above and below by the lines between points 1 and 3 and 2 and 3,
      # and horizontally between the vertical left side formed by points 1 and 2, and point 3 on the right
      # where the upper and lower bounding lines meet.
      const slope_23 = slope_between(pts[2], pts[3])
      integral_fn_between_lines_meeting_at_point_and_vert_line(f_x_mon, pts[3], slope_13, slope_23, pts[1]._1, mesh)
    else
      # Points 1 and 2 do not form a vertical line. Integrate between points 1 and 2, and between points 2 and 3 if
      # points 2 and 3 don't lie on a vertical line.
      const fst_seg = integral_fn_between_lines_meeting_at_point_and_vert_line(f_x_mon, pts[1], slope_12, slope_13, pts[2]._1, mesh)
      const snd_seg =
        if pts[2]._1 == pts[3]._1 # vertical right side, no second segment
          zeroR
        else
          const slope_23 = slope_between(pts[3], pts[2])
          integral_fn_between_lines_meeting_at_point_and_vert_line(f_x_mon, pts[3], slope_13, slope_23, pts[2]._1, mesh)
        end
      fst_seg + snd_seg
    end
  else # side face
    # We want to compute int_{t=0..1} f(p(t)) mon(p(t)-o) |p'(t)| dt, where p:[0,1] -> R^2
    # is a bijection traversing the side face smoothly, p(0) and p(1) being the side
    # endpoints, and o is the local origin for the side, which is the side's midpoint.
    const a,b = side_face_endpoint_pair(face,
                                       fe_v1, fe_v1 + ref_tri.v12, fe_v1 + ref_tri.v13,
                                       ref_tri.nums_side_faces_between_vertexes,
                                       false) # lesser endpoints first => false
    const o = 1/2*(a + b) # The local origin of each side face is the midpoint.
    const side_len = hypot(b._1 - a._1, b._2 - a._2)

    const integrand = let p_t = mesh.integrand_work_array
      function(t_::Vector{R}) # compute f(p(t)) mon(p(t)-o) |p'(t)|
        const t = t_[1]
        p_t[1] = a._1 + (b._1 - a._1)*t
        p_t[2] = a._2 + (b._2 - a._2)*t
        const f_val = f(p_t)
        # Make p_t side origin relative for the monomial evaluation.
        p_t[1] -= o._1; p_t[2] -= o._2
        const mon_val = Poly.value_at(mon, p_t)
        f_val * mon_val * side_len  # |p'(t)| = |b-a| = side length
      end
    end

    hcubature(integrand,
              singleZero, singleOne,
              reltol=mesh.integration_rel_err, abstol=mesh.integration_abs_err)[1]
  end
end

# Integrate a side-local monomial vs. a vector monomial interpreted relative to the
# entire finite element, dot multiplied with the outward normal for the side.
import Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side
function integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(mon::Monomial,
                                                                     vmon::VectorMonomial,
                                                                     oshape::OShape, side_face::FEFaceNum,
                                                                     mesh::TriMesh)
  const ref_tri = mesh.oshapes[oshape]

  # fe-relative side endpoints
  const a,b = side_face_endpoint_pair(side_face,
                                     zeroPt, ref_tri.v12, ref_tri.v13,
                                     ref_tri.nums_side_faces_between_vertexes,
                                     false) # lesser endpoints first => false
  # The local origin of each side face is the midpoint.
  const side_o = 1/2*(a + b)
  const onormal = ref_tri.outward_normals_by_side_face[side_face]
  const side_len = hypot(b._1 - a._1, b._2 - a._2)

  # For path p:[0,1] -> S traversing from a to b linearly in fe-relative coordinates, compute
  #   mon(p(t)-side_o) vmon(p(t)).n  |p'(t)|
  const integrand = let p_t = mesh.integrand_work_array
    function(t_::Vector{R})
      const t = t_[1]
      # fe-relative path point
      p_t[1] = a._1 + (b._1 - a._1)*t
      p_t[2] = a._2 + (b._2 - a._2)*t
      const vmon_dot_normal = let vmon_mon_val = Poly.value_at(vmon.mon, p_t)
        vmon.mon_pos == 1 ? onormal._1 * vmon_mon_val: onormal._2 * vmon_mon_val
      end
      # Make p_t relative to the side origin for the monomial evaluation.
      p_t[1] -= side_o._1; p_t[2] -= side_o._2
      const mon_val = Poly.value_at(mon, p_t)
      mon_val * vmon_dot_normal * side_len
    end
  end

  # TODO: really should do this integral analytically
  hcubature(integrand,
            singleZero, singleOne,
            reltol=mesh.integration_rel_err, abstol=mesh.integration_abs_err)[1]
end

# Integrate a finite element relative monomial vs. a side relative monomial.
import Mesh.integral_fe_rel_x_side_rel_on_oshape_side
function integral_fe_rel_x_side_rel_on_oshape_side(fe_rel_mon::Monomial,
                                                   side_rel_mon::Monomial,
                                                   oshape::OShape, side_face::FEFaceNum,
                                                   mesh::TriMesh)
  const ref_tri = mesh.oshapes[oshape]

  # fe-relative side endpoints
  const a,b = side_face_endpoint_pair(side_face,
                                     zeroPt, ref_tri.v12, ref_tri.v13,
                                     ref_tri.nums_side_faces_between_vertexes,
                                     false) # lesser endpoints first => false
  # The local origin of each side face is the midpoint.
  const side_o = 1/2*(a + b)
  const onormal = ref_tri.outward_normals_by_side_face[side_face]
  const side_len = hypot(b._1 - a._1, b._2 - a._2)

  # For path p:[0,1] -> S traversing from a to b linearly in fe-relative coordinates, compute
  #   side_rel_mon(p(t)-side_o) fe_rel_mon(p(t))  |p'(t)|
  const integrand = let p_t = mesh.integrand_work_array
    function(t_::Vector{R})
      const t = t_[1]
      # fe-relative path point
      p_t[1] = a._1 + (b._1 - a._1)*t
      p_t[2] = a._2 + (b._2 - a._2)*t
      const fe_rel_mon_val = Poly.value_at(fe_rel_mon, p_t)
      # Make p_t relative to the side origin for the side-local monomial evaluation.
      p_t[1] -= side_o._1; p_t[2] -= side_o._2
      const side_rel_mon_val = Poly.value_at(side_rel_mon, p_t)
      fe_rel_mon_val * side_rel_mon_val * side_len
    end
  end

  # TODO: really should do this integral analytically
  hcubature(integrand,
            singleZero, singleOne,
            reltol=mesh.integration_rel_err, abstol=mesh.integration_abs_err)[1]
end

# required AbstractMesh functions
####################################################################


# Integration Auxiliary Functions

# Integrate a monomial over the triangular region bounded on two sides by two
# non-vertical lines of indicated slopes which meet at a point p, and by the
# indicated vertical line as the remaining side.
function integral_mon_between_lines_meeting_at_point_and_vert_line(mon::Monomial,
                                                                   p::Point,
                                                                   slope_1::R, slope_2::R,
                                                                   vert_line_x::R)
  # slopes of the lower and upper lines
  const m_l, m_u = p._1 < vert_line_x ? (min(slope_1, slope_2), max(slope_1, slope_2)) :
                                        (max(slope_1, slope_2), min(slope_1, slope_2))
  # polynomials representing the lower and upper lines for the inner (y) integral's bounds
  const y_l, y_u = let x = MON_VAR_1D
    p._2 + m_l*(x-p._1),
    p._2 + m_u*(x-p._1)
  end

  # Compute inner integral, which is an integral over y with x held constant, as a polynomial in x.
  const mon_antider = Poly.antideriv(dim(2), mon)
  const inner_int = Poly.reduce_dim_by_subst(dim(2), y_u, mon_antider) - Poly.reduce_dim_by_subst(dim(2), y_l, mon_antider)

  # Difference of anti-derivative of inner integral over 1st coordinate bounds is the final integral value.
  const antider_inner_int = Poly.antideriv(dim(1), inner_int)
  # bounds of integration for outer integral
  const ib_l, ib_u = p._1 < vert_line_x ? (p._1, vert_line_x) : (vert_line_x, p._1)
  Poly.value_at(antider_inner_int, ib_u) - Poly.value_at(antider_inner_int, ib_l)
end

# Integrate a function over the triangular region bounded on two sides by two
# non-vertical lines of indicated slopes which meet at a point p, and by the
# indicated vertical line as the remaining side.
function integral_fn_between_lines_meeting_at_point_and_vert_line(g::Function,
                                                                  p::Point,
                                                                  slope_1::R, slope_2::R,
                                                                  vert_line_x::R,
                                                                  mesh::TriMesh)
  const xmin, xmax = min(p._1, vert_line_x), max(p._1, vert_line_x)
  const w = xmax - xmin # width of triangular integration region

  # Method
  # We will pull back the integration over the original triangular section T to an integration
  # over the unit square, by change of variables via the bijection
  #   t:[0,1]^2 -> T
  #   t(x,y) = (xmin + x w,  p_2 + m1 (xmin + x w - p_1) + y(m2 - m1)(xmin + x w - p_1))
  #     for (x,y) in [0,1]^2.  Here m1 and m2 are the slopes of the non-vertical bounding lines.
  # The determinant of the derivative matrix Dt(x,y) is
  #                     |          w                           0               |
  #  det Dt(x,y)| = det |                                                      |
  #                     | m1 w + y(m2 - m1) w      (m2 - m1)(xmin + x w - p_1) |
  #               = w (m2 - m1) (xmin + x w - p_1).
  # Now by applying change of variables in the integral of g over the T via the mapping t,
  # we have
  #   int_T g = int_0^1 int_0^1 g(t(x,y)) |det Dt(x,y)| dy dx
  #            = int_0^1 int_0^1 g(t(x,y)) |w (m2 - m1) (xmin + x w - p_1)| dy dx

  const slopediff = slope_2 - slope_1
  const w_slopediff = w * slopediff

  function gt_absdetDt(r::Vector{R}) # evaluate integrand at rectangle cubature point r
    const x, y = r[1], r[2]
    const xT = xmin + x * w      # (triangle x)
    const xT_minus_p1 = xT - p._1
    try
      # Construct the triangle point t(x,y).
      const t = r # reuse the rectangle point array temporarily for the triangle point
      t[1] = xT
      t[2] = p._2 + slope_1 * xT_minus_p1  +  y * slopediff * xT_minus_p1
      const det_Dt = w_slopediff * xT_minus_p1
      g(t) * abs(det_Dt)
    finally
      # restore input array to its original state
      r[1] = x; r[2] = y
    end
  end

  hcubature(gt_absdetDt,
            mesh.space_dim_zeros, mesh.space_dim_ones, # unit square bounds
            reltol=mesh.integration_rel_err, abstol=mesh.integration_abs_err)[1]
end


###############
# Etc

const SPACE_DIM= dim(2)
const zeroPt = Point(0.,0.)
const singleZero = [zeroR]
const singleOne = [oneR]
const MON_VAR_1D = Monomial(1)

# Reference triangle for a finite element.
ref_tri_for_fe(fe::FENum, mesh::TriMesh) = mesh.oshapes[mesh.fes[fe].oshape]

physical_region_tag(fe::FENum, mesh::TriMesh) = mesh.phys_reg_tags_by_FENum[fe]

geometric_entity_tag(fe::FENum, mesh::TriMesh) = mesh.geom_ent_tags_by_FENum[fe]


# vector operations

slope_between(p::Point, q::Point) = (q._2 - p._2)/(q._1 - p._1)

*/

enum SideEndpointsOrdering {
  LesserEndpointsFirst,
  FirstTraversedEndpointsFirst,
}

fn side_face_endpoint_pair(sf: SideFace,
                           v0: Point, v1: Point, v2: Point,
                           nums_side_faces_btw_verts: (u8,u8,u8),
                           endpoints_ordering: SideEndpointsOrdering) -> (Point, Point)
{
  let mut cur_sf = 0u; 
  
  let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
  let side_faces = [(v0, v1, sfs_v01),
                    (v1, v2, sfs_v12),
                    (v2, v0, sfs_v20)];
  
  for &(va, vb, num_sfs_btw_va_vb) in side_faces.iter()
  {
    match num_sfs_btw_va_vb
    {
      1 => {
        if cur_sf == *sf { return mk_endpoint_pair(va, vb, endpoints_ordering) }
        cur_sf += 1;
      }
      2 => {
        let midpt = { let ((va_x,va_y),(vb_x,vb_y)) = (va,vb); (0.5*(va_x + vb_x), 0.5*(va_y + vb_y)) };
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

fn mk_endpoint_pair(pt1: Point, pt2: Point, endpoints_ordering: SideEndpointsOrdering) -> (Point, Point)
{
  match endpoints_ordering
  {
    LesserEndpointsFirst => if pt1 < pt2 { (pt1, pt2) } else { (pt2, pt1) },
    FirstTraversedEndpointsFirst => (pt1, pt2)
  }
}



fn norm((v_x, v_y): Vec) -> R { hypot(v_x, v_y) }


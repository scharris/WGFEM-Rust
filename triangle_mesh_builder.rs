use triangle_mesh::*;
use common::*;
use monomial::{Monomial};
use mesh::*;
use storage_by_ints::StorageByInts2;

use std::vec;
use std::num::{pow_with_uint};
use std::iter::{Iterator, range_inclusive};
use std::hashmap::{HashMap, HashSet};
use std::util::swap;

// Mesh Builder
pub struct MeshBuilder {
  
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
  
  pub fn make_mesh<Mon:Monomial,I:Iterator<BaseTri>>(base_pts_by_nodenum: ~[Point],
                                                     mut base_tris_iter: &mut I,
                                                     est_base_tris: uint,
                                                     global_subdiv_iters: uint,
                                                     integration_rel_err: R,
                                                     integration_abs_err: R,
                                                     load_phys_reg_tags: bool,
                                                     load_geom_ent_tags: bool) -> TriMesh<Mon> {

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
    for base_tri in base_tris_iter {
      bldr.process_base_tri(&base_tri, base_pts_by_nodenum, global_subdiv_iters);
    }

    // Create the final non-boundary sides data structures based on the mapping of side endpoints to fe faces.
    let (nbsidenums_by_fe_face, nbsideincls_by_nbsidenum) = bldr.create_nb_sides_data();
    assert!(nbsideincls_by_nbsidenum.len() == bldr.num_nb_sides as uint);

    TriMesh::new(bldr.take_fes(),
                 bldr.take_oshapes(),
                 nbsidenums_by_fe_face,
                 nbsideincls_by_nbsidenum,
                 bldr.num_b_sides as uint,
                 integration_rel_err,
                 integration_abs_err,
                 bldr.take_phys_reg_tags(),
                 bldr.take_geom_ent_tags())
  }

  fn process_base_tri(& mut self,
                      base_tri: &BaseTri,
                      mesh_pts_by_nodenum: &[Point],
                      global_subdiv_iters: uint) {

    let (v0, v1, v2) = (mesh_pts_by_nodenum[*base_tri.node_num_0],
                        mesh_pts_by_nodenum[*base_tri.node_num_1],
                        mesh_pts_by_nodenum[*base_tri.node_num_2]);

    // For each of the three sides of this mesh element to be subdivided, the generated subdivision elements
    // can be made to have more than one face between their vertexes lying on the side.  This is useful to
    // support "hanging" nodes where a finer subdivision is adjacent to this one.  The triplet returned
    // represents the number of side faces that a generated element has if its vertexes lie between base
    // triangle vertex pairs v0 and v1, v1 and v2, and v2 and v0, respectively.
    let nums_sfs_btw_verts = base_tri.nums_side_faces_between_vertexes();
    
    // Add any extra subdivision iterations to be done in this element.
    let subdiv_iters = global_subdiv_iters + base_tri.extra_subdiv_iters();

    // Register the primary reference triangles for our mesh element's subdivisions. Normally there will only be one
    // such primary (ie. non-inverted) reference triangle, however when multiple side faces are present between two
    // vertexes then additional reference triangles are required.
    let first_new_pri_oshape = self.register_primary_ref_tris(v0,v1,v2, nums_sfs_btw_verts, subdiv_iters);
    let last_new_pri_oshape = OShape(self.oshapes.len()-1);

    // Create a lookup function for the new primary oshape numbers by their numbers of side faces between vertexes.
    let pri_oshapes_by_nums_sfs_btw_verts: |(u8,u8,u8)| -> OShape = |sfs_btw_verts: (u8,u8,u8)| {
      for os in range_inclusive(*first_new_pri_oshape, *last_new_pri_oshape) {
        if self.oshapes[os].nums_side_faces_between_vertexes == sfs_btw_verts {
          return OShape(os);
        }
      }
      fail!("Reference triangle not found by numbers of side faces between vertexes.");
    };

    if subdiv_iters > 0 {
      // We're subdividing, so register the secondary reference triangle.
      let sec_oshape = self.register_secondary_ref_tri(v0,v1,v2, subdiv_iters);

      // Do the subdivisions.
      self.subdivide_primary(v0, v1, v2,
                             subdiv_iters,
                             nums_sfs_btw_verts,
                             pri_oshapes_by_nums_sfs_btw_verts, sec_oshape,
                             base_tri.tag_phys_reg, base_tri.tag_geom_ent);
    } else { // no subdivision to be done
      // The base triangle itself is our finite element, with its own reference triangle.
      self.add_fe(first_new_pri_oshape, v0,v1,v2, base_tri.tag_phys_reg, base_tri.tag_geom_ent);
    }
  }
 
 /* Create primary reference triangles for the given triangle to be subdivided. These reference triangles
  * are rescaled translations of the original undivided triangle. If nums_side_faces_btw_verts is other
  *  than (1,1,1), then multiple reference triangles will be generated, differing in the number of side
  *  faces between their vertexes for support of "hanging" nodes.  The function returns a function mapping
  *  the numbers of side faces between vertexes to the newly registered reference triangle (oshape) number.
  *  Pains are taken to not create more reference triangles than are actually used, so that some global
  *  mesh properties such as maximum element diameter can be determined from only the reference elements.
  *  Returns the first new oriented shape number that was created.
  */
  fn register_primary_ref_tris(& mut self,
                               v0: Point, v1: Point, v2: Point,
                               nums_sfs_btw_verts: (u8,u8,u8),
                               subdiv_iters: uint) -> OShape {

    let first_new_oshape = OShape(self.oshapes.len());
   
    let scale = { let p: uint = pow_with_uint(2u, subdiv_iters); 1./(p as R) };

    if subdiv_iters == 0 || nums_sfs_btw_verts == (1u8, 1u8, 1u8) {
      self.oshapes.push(RefTri::new(v0,v1,v2, scale, nums_sfs_btw_verts));
    } else { // We are subdividing, and there is more than one side face between some vertex pair.
      let (nums_sfs_btw_verts_0, nums_sfs_btw_verts_1, nums_sfs_btw_verts_2) = nums_sfs_btw_verts;

      // We need to add a separate primary reference triangle for each triplet representing the number of side
      // faces between vertexes that occur in primary subdivision triangles of this base triangle.  We use the
      // set below to track the unique triplets of numbers of side faces between vertexes.
      self.sfs_btw_verts_work_set.clear();
    
      // Handle subdivision triangles at one of the base triangle's corners, v0, v1 or v2.
      // To support the corner elements, we need primary reference triangles which get 2 of their 3 numbers
      // of faces between vertex pairs from nums_side_faces_btw_verts, because they have two sides on the
      // original undivided triangle's sides.
      self.sfs_btw_verts_work_set.insert((nums_sfs_btw_verts_0, 1u8, nums_sfs_btw_verts_2)); // v0 corner element
      self.sfs_btw_verts_work_set.insert((nums_sfs_btw_verts_0, nums_sfs_btw_verts_1, 1u8)); // v1 corner element
      self.sfs_btw_verts_work_set.insert((1u8, nums_sfs_btw_verts_1, nums_sfs_btw_verts_2)); // v2 corner element

      // If subdividing more than once, then we also need:
      //  - a primary reference triangle with only one side face between each of its vertex pairs (within the central
      //    secondary triangle).
      //  - primary reference triangles which only have one pair of vertexes contained in one of the original base
      //    triangle's sides, so that one pair of vertexes may have multiple side faces between them but not so for
      //    the other two vertex pairs.
      if subdiv_iters > 1 {
        self.sfs_btw_verts_work_set.insert((1u8, 1u8, 1u8));
        self.sfs_btw_verts_work_set.insert((nums_sfs_btw_verts_0, 1u8, 1u8));
        self.sfs_btw_verts_work_set.insert((1u8, nums_sfs_btw_verts_1, 1u8));
        self.sfs_btw_verts_work_set.insert((1u8, 1u8, nums_sfs_btw_verts_2));
      } 

      // Now create the reference triangles for each unique triplet of numbers of side faces between vertexes.
      for &sfs_btw_verts in self.sfs_btw_verts_work_set.iter() {
        self.oshapes.push(RefTri::new(v0,v1,v2, scale, sfs_btw_verts));
      }
    }
    
    first_new_oshape 
  }

  // Register the single secondary (inverted) subdivision triangle for the given (primary) base triangle, and number of
  // subdivisions to be applied on the base (primary) triangle.
  fn register_secondary_ref_tri(& mut self,
                                pv0: Point, pv1: Point, pv2: Point, // base (primary) triangle vertexes
                                subdiv_iters: uint) -> OShape {
    let sec_oshape = OShape(self.oshapes.len());
    let sec_ref_tri = {
      let scale = { 
        let p: uint = pow_with_uint(2u, subdiv_iters-1); // midpoints already represent one subdivision
        1./(p as R)
      }; 
      RefTri::new(midpt(pv0, pv1), midpt(pv1, pv2), midpt(pv2, pv0), scale, (1u8,1u8,1u8))
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
                       tag_phys_reg: Tag, tag_geom_ent: Tag) {
    if iters == 0 {
      let oshape = pri_oshapes_by_nums_sfs_btw_verts(nums_side_faces_btw_verts);
      self.add_fe(oshape, v0,v1,v2, tag_phys_reg, tag_geom_ent);
    } else {
      let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
      let (midpt_v01, midpt_v12, midpt_v20) = (midpt(v0,v1), midpt(v1,v2), midpt(v2,v0));

      // The sub-triangle including v0 has its first and last vertex pairs embedded in the original triangle's first
      // and third sides, and so inherits the corresponding numbers of faces between vertexes from the first and third
      // entries of nums_side_faces_btw_verts. The middle vertex pair (midpt_v01 to midpt_v20) of this sub-triangle
      // does not lie along an original side, however, so has only one face between these vertexes.  Similar arguments
      // apply for the remaining sub-triangles.
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
                         tag_phys_reg: Tag, tag_geom_ent: Tag) {
    if iters == 0 { 
      self.add_fe(sec_oshape, v0,v1,v2, tag_phys_reg, tag_geom_ent);
    } else {
      let (midpt_v01, midpt_v12, midpt_v20) = (midpt(v0,v1), midpt(v1,v2), midpt(v2,v0));

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

      // The central sub-triangle is a primary triangle.  We need to order the vertexes properly, so that all primary
      // and secondary triangles will have the same vertex orientations, in the sense that for any two subtriangles of
      // the same type (primary or secondary), some translation can take the vertexes of one to the corresponding
      // vertexes of the other. We do this by numbering the primary vertexes just as they were numbered in the original
      // undivided primary triangle, and numbering all secondary vertexes just as in the original central secondary
      // triangle of the undivided primary triangle.
      // We shifted "forward" when numbering vertexes of a central secondary subtriangle from the enclosing primary
      // triangle, so that e.g. the secondary's v0 is the midpoint of enclosing primary's v0 and v1. So now we must
      // shift vertex numbers "back" for this enclosed primary triangle, so e.g. the central primary triangle's v0
      // is the midpoint of this enclosing secondary triangle's v2 and v0, and similarly for the other vertexes.
      self.subdivide_primary(midpt_v20, midpt_v01, midpt_v12,
                             iters-1,
                             (1u8,1u8,1u8),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag_phys_reg, tag_geom_ent);
    }
  }

  // This procedure will be used to create all finite elements during the processing of the input base triangles.
  fn add_fe(& mut self, oshape: OShape, v0: Point, v1: Point, v2: Point, tag_phys_reg: Tag, tag_geom_ent: Tag) {

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
                                        nums_side_faces_btw_verts: (u8,u8,u8)) -> (int, int) {
    let mut nb_sides_delta = 0;
    let mut b_sides_delta = 0;
    
    // Take the endpoint pairs work buffer, filled with the side face endpoint pairs.
    let endpt_pairs = MeshBuilder::fill_side_face_endpoint_pairs(self.take_side_face_endpoint_pairs_buf(),
                                                                 v0,v1,v2, 
                                                                 nums_side_faces_btw_verts,
                                                                 LesserEndpointsFirst);

    // Register the side face representations for this element by their endpoints.
    for sf in range(0, endpt_pairs.len()) {
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

 /* This function defines the side faces enumeration for a finite element of given vertexes and numbers of faces
  * between vertexes. Side faces are returned as an array of side endpoint pairs indexed by side face number. If
  * lesser_endpts_first is true, then each endpoint pair endpoint will have the lesser point (compared
  * lexicographically) in the first component of the pair.
  */
  fn fill_side_face_endpoint_pairs(mut endpts_buf: ~[(Point,Point)],
                                   v0: Point, v1: Point, v2: Point,
                                   nums_side_faces_btw_verts: (u8,u8,u8),
                                   endpoints_ordering: SideEndpointsOrdering) -> ~[(Point,Point)] {
    endpts_buf.clear();

    let (sfs_v01, sfs_v12, sfs_v20) = nums_side_faces_btw_verts;
    let side_faces = [(v0, v1, sfs_v01),
                      (v1, v2, sfs_v12),
                      (v2, v0, sfs_v20)];

    for &(va, vb, num_sfs_btw_va_vb) in side_faces.iter() {
      match num_sfs_btw_va_vb {
        1 => endpts_buf.push(mk_endpoint_pair(va, vb, endpoints_ordering)),
        2 => {
          let midpt = midpt(va, vb);
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
  fn create_nb_sides_data(& mut self) -> (StorageByInts2<Option<NBSideNum>>, ~[NBSideInclusions]) {

    let max_sides_per_fe = self.oshapes.iter().map(|rt| rt.num_side_faces).max().unwrap();
    
    let mut nbsidenums_by_fe_face = StorageByInts2::from_elem(self.fes.len(), max_sides_per_fe, None);
    let mut nbsideincls_by_nbsidenum = vec::with_capacity(self.num_nb_sides);

    for (_, side_reps) in self.side_reps_by_endpts.iter() {
      match side_reps.snd_rep {
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

  fn take_fes(& mut self) -> ~[ElTri] {
    let mut fes = ~[];
    swap(& mut fes, & mut self.fes);
    fes
  }
  
  fn take_oshapes(& mut self) -> ~[RefTri] {
    let mut oshapes = ~[];
    swap(& mut oshapes, & mut self.oshapes);
    oshapes
  }
  
  fn take_phys_reg_tags(& mut self) -> Option<~[Tag]> {
    let mut tags = None;
    swap(& mut tags, & mut self.phys_reg_tags_by_fenum);
    tags 
  }
  
  fn take_geom_ent_tags(& mut self) -> Option<~[Tag]> {
    let mut tags = None;
    swap(& mut tags, & mut self.geom_ent_tags_by_fenum);
    tags 
  }

  fn take_side_face_endpoint_pairs_buf(& mut self) -> ~[(Point,Point)] {
    let mut endpt_pairs_opt = None;
    swap(& mut endpt_pairs_opt, & mut self.sf_endpt_pairs_work_buf);
    endpt_pairs_opt.unwrap()
  }
  
  fn return_side_face_endpoint_pairs_buf(& mut self, buf: ~[(Point,Point)]) {
    let mut ret_buf_opt = Some(buf);
    swap(& mut ret_buf_opt, & mut self.sf_endpt_pairs_work_buf);
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
  */

} // MeshBuilder impl


// A structure to hold the one or two local representations of any side as fe/side face pairs.
struct SideReps {
  fst_rep: (FENum, SideFace),
  snd_rep: Option<(FENum, SideFace)>,
}


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
  
  fn extra_subdiv_iters(& self) -> uint {
    0u // TODO: Allow specifying extra iterations via a Gmsh tag.
  }

 /* For each of the three sides of the indicated base mesh triangle to be subdivided, we can specify here the number
  * of faces that the generated subtriangles should have between their vertexes which lie on this side. This allows
  * these subdivision elements to meet those of a finer subdivision in a mesh element adjacent to this element
  * ("hanging nodes").  The numbers are returned as a triplet of integers, corresponding to the face counts to be
  * generated for elements' sides within sides v0v1, v1v2, and v2v0.
  */
  fn nums_side_faces_between_vertexes(& self) -> (u8,u8,u8) {
    (1u8,1u8,1u8) // TODO: Allow specifying somewhere such as in mesh element tags.
  }
  
} // BaseTri impl

impl SideReps {
  
  fn new(fe_sf: (FENum,SideFace)) -> SideReps { SideReps { fst_rep: fe_sf, snd_rep: None } } 

  fn add(& mut self, (fe,sf): (FENum,SideFace)) {
    // Keep the representations in fe order.
    match self.fst_rep {
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
# Construct from Gmsh .msh formatted input stream.
function TriMesh(is::IO, subdiv_iters::Integer,
                 integration_rel_err::R = 1e-12, integration_abs_err::R = 1e-12,
                 load_phys_reg_tags::Bool = false, load_geom_ent_tags::Bool = false)

  const base_pts_by_nodenum = read_gmsh_nodes(is)

  # Read to the beginning of the elements section.
  if !read_through_line_starting(is, "\$Elements")
    error("Could not find beginning of elements section in mesh file.")
  end

  # Estimate number of input mesh triangles from the number of elements declared and the number
  # of skipped lower order elements (points and lines) before the first triangle element.
  const decl_num_els = uint64(readline(is)) # This count can include unwanted lower order elements.
  const line, lines_read = read_until(is, is_polytope_or_endmarker)
  if line == "" || beginswith(line, "\$EndElements")
    error("No elements found in \$Elements section.")
  end
  const est_base_tris = decl_num_els - (lines_read - 1)

  const base_tris_iter = GmshElementsIter(line, is)

  TriMesh(base_pts_by_nodenum,
          base_tris_iter,
          est_base_tris,
          subdiv_iters,
          integration_rel_err, integration_abs_err,
          load_phys_reg_tags, load_geom_ent_tags)
end





# ----------------------------------
# Gmsh msh format reading support

typealias ElTypeNum Uint8
eltypenum(i::Integer) = if i>=1 && i<=typemax(ElTypeNum) convert(ElTypeNum, i) else error("invalid element type: $i") end
eltypenum(s::String) = eltypenum(int(s))

# Gmsh triangle type codes
const ELTYPE_3_NODE_TRIANGLE = eltypenum(2)
# lower order elements, which can be ignored
const ELTYPE_POINT = eltypenum(15)
const ELTYPE_2_NODE_LINE = eltypenum(1)
const ELTYPE_3_NODE_LINE = eltypenum(8)
const ELTYPE_4_NODE_LINE = eltypenum(26)
const ELTYPE_5_NODE_LINE = eltypenum(27)
const ELTYPE_6_NODE_LINE = eltypenum(28)

function is_3_node_triangle_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_3_NODE_TRIANGLE
end

function is_lower_order_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_POINT ||
  el_type == ELTYPE_2_NODE_LINE ||
  el_type == ELTYPE_3_NODE_LINE ||
  el_type == ELTYPE_4_NODE_LINE ||
  el_type == ELTYPE_5_NODE_LINE ||
  el_type == ELTYPE_6_NODE_LINE
end

const TOKN_NODELINE_ELNUM = 1
const TOKN_NODELINE_POINT1 = 2
const TOKN_NODELINE_POINT2 = 3
const TOKN_NODELINE_POINT3 = 4

# Gmsh element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELLINE_ELTYPE = 2
const TOKN_ELLINE_NUMTAGS = 3


# Gmsh base elements iterator

immutable GmshElementsIter
  first_base_tri_line::String
  is::IO # stream should be positioned after first base element line
end

import Base.start, Base.done, Base.next
function start(gmsh_els_iter::GmshElementsIter)
  base_tri_from_gmsh_el(split(gmsh_els_iter.first_base_tri_line, ' '))
end

function done(gmsh_els_iter::GmshElementsIter, next_base_tri)
  next_base_tri == nothing
end

function next(gmsh_els_iter::GmshElementsIter, next_base_tri)
  while true
    const line = readline(gmsh_els_iter.is)
    if beginswith(line, "\$EndElements")
      return (next_base_tri, nothing)
    elseif line == ""
      error("No EndElements marker found.")
    else
      const toks = split(line, ' ')
      const el_type = eltypenum(toks[TOKN_ELLINE_ELTYPE])
      if is_lower_order_el_type(el_type)
        continue
      elseif is_3_node_triangle_el_type(el_type)
        return (next_base_tri, base_tri_from_gmsh_el(toks))
      else
        error("Element type $el_type is not supported.")
      end
    end
  end
end

const EMPTY_OTHER_TAGS = Array(Tag,0)

function base_tri_from_gmsh_el(toks::Array{String,1})
  const num_tags = int(toks[TOKN_ELLINE_NUMTAGS])
  assert(num_tags >= 2)
  const physreg_tag = tag(int64(toks[TOKN_ELLINE_NUMTAGS+1]))
  const geoment_tag = tag(int64(toks[TOKN_ELLINE_NUMTAGS+2]))
  const other_tags = num_tags >= 3 ? map(t -> tag(int64(t)), toks[TOKN_ELLINE_NUMTAGS+3:TOKN_ELLINE_NUMTAGS+num_tags]) :
                                     EMPTY_OTHER_TAGS
  const point_nums = let last_tag_tokn = TOKN_ELLINE_NUMTAGS + num_tags;
    (nodenum(toks[last_tag_tokn + 1]),
     nodenum(toks[last_tag_tokn + 2]),
     nodenum(toks[last_tag_tokn + 3]))
  end

  BaseTri(point_nums,
          physreg_tag,
          geoment_tag,
          other_tags)
end


###############################################
# Mesh File Reading Utilities

function read_gmsh_nodes(is::IO)
  if !read_through_line_starting(is, "\$Nodes")
    error("Nodes section not found in mesh input file.")
  else
    # Next line should be the vert count.
    const count = int(readline(is))
    const pts = Array(Point, count)
    l = readline(is)
    while !beginswith(l, "\$EndNodes") && l != ""
      const toks = split(l, ' ')
      const pt = Point(convert(R, float64(toks[TOKN_NODELINE_POINT1])), convert(R, float64(toks[TOKN_NODELINE_POINT2])))
      if length(toks) >= TOKN_NODELINE_POINT3 && float64(toks[TOKN_NODELINE_POINT3]) != 0.0
        error("Nodes with non-zero third coordinates are not supported in this 2D mesh reader.")
      end
      pts[uint64(toks[TOKN_NODELINE_ELNUM])] = pt
      l = readline(is)
    end
    if !beginswith(l, "\$EndNodes") error("End of nodes section not found in mesh file.") end
    pts
  end
end

function is_polytope_or_endmarker(line::ASCIIString)
  if beginswith(line, "\$EndElements")
    true
  else
    const toks = split(line, ' ')
    !is_lower_order_el_type(eltypenum(toks[TOKN_ELLINE_ELTYPE]))
  end
end

function read_through_line_starting(io::IO, line_start::ASCIIString)
  l = readline(io)
  while l != "" && !beginswith(l, line_start)
    l = readline(io)
  end
  l != ""
end

# Read the stream until the line condition function returns true or the stream is exhausted, returning
# either the line passing the condition, or "" if the stream was exhausted, and (in either case)
# the number of lines read including the successful line if any.
function read_until(io::IO, line_cond_fn::Function)
  l = readline(io)
  lines_read = 1
  while l != "" && !line_cond_fn(l)
    l = readline(io)
    lines_read += 1
  end
  l, lines_read
end

# Mesh File Reading Utilities
###############################################



# Gmsh msh format reading support
# ----------------------------------
*/

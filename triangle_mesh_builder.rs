use triangle_mesh::*;
use common::*;
use monomial::{Monomial};
use mesh::*;
use storage_by_ints::StorageByInts2;

use std::vec;
use std::num::{pow_with_uint, from_uint};
use std::iter::{Iterator, range_inclusive};
use std::hashmap::{HashMap, HashSet};
use std::util::swap;
use std::io::buffered::BufferedReader;


// Triangle mesh builder type.
pub struct TriMeshBuilder {
  
  oshapes: ~[RefTri],
  fes: ~[ElTri],

  side_reps_by_endpts: HashMap<(Point,Point),SideReps>,

  tags_by_fenum: Option<~[Tag]>,

  num_nb_sides: uint,
  num_b_sides: uint,

  // work buffers
  sfs_btw_verts_work_set: HashSet<(u8,u8,u8)>,
  sf_endpt_pairs_work_buf: Option<~[(Point,Point)]>,

}


impl TriMeshBuilder {

  // Construct a triangle mesh from a Gmsh msh formatted stream.
  pub fn from_gmsh_msh_stream <Mon: Monomial, I: Reader>
         ( msh_is: &mut BufferedReader<I>,
           subdiv_iters: uint,
           integration_rel_err: R, integration_abs_err: R,
           load_tags: bool )
         -> TriMesh<Mon>
  {
    let base_pts_by_nodeix = read_msh_nodes(msh_is);

    // Iterate lines beginning at the Elements section.
    let mut lines = msh_is.lines().skip_while(|line| !line.starts_with("$Elements"));

    if lines.next().is_none() {
      fail!("$Elements section not found in msh file.");
    }

    // Estimate number of input mesh triangles from the number of elements declared, then read to the first triangle element,
    // skipping any leading lower order elements and subtracting their count from the estimated number of elements.
    let num_els_line = lines.next().unwrap();
    let mut est_base_tris: uint = from_str(num_els_line.trim_right()).unwrap(); // This count can include unwanted lower order elements.
    let mut first_tri_line = None;
    lines.advance(|line| {
      let el_type: uint = extract_field(line, TOKIX_ELLINE_ELTYPE);
      if is_lower_order_el_type(el_type) {
        if est_base_tris > 0 { est_base_tris -= 1; }
        true
      } else {
        first_tri_line = Some(line);
        false
      }
    });

    if first_tri_line.is_none() {
      fail!("No triangle elements found.");
    }
 
    // Read the base triangles.
    let mut base_tris = vec::with_capacity(est_base_tris);
    for line in first_tri_line.move_iter().chain(lines.take_while(|line| !line.starts_with("$EndElements"))) {
      match base_tri_from_line(line, base_pts_by_nodeix) {
        Some(base_tri) => base_tris.push(base_tri),
        _ => {}
      }
    }

    TriMeshBuilder::from_base_tris(base_tris,
                                   subdiv_iters,
                                   integration_rel_err, integration_abs_err,
                                   load_tags)
  }


  pub fn from_base_tris<Mon: Monomial>
         ( base_tris: &[BaseTri],
           global_subdiv_iters: uint,
           integration_rel_err: R, integration_abs_err: R,
           load_tags: bool )
         -> TriMesh<Mon>
  {
    // Make estimates of the number of finite elements and reference triangles for storage allocation.
    // This estimate ignores optional additional subdivisions that may be specified for some input elements.
    let est_fes = base_tris.len() * pow_with_uint(4, global_subdiv_iters);
    // We'll estimate 2 reference triangles for each mesh element if we are subdividing, 1 if not.
    // We ignore for this estimate any additional reference triangles for additional subdivisions.
    let est_ref_tris = base_tris.len() * (if global_subdiv_iters > 0 { 2 } else { 1 });
    
    let mut bldr = TriMeshBuilder {
      fes: vec::with_capacity(est_fes),
      oshapes: vec::with_capacity(est_ref_tris),
      side_reps_by_endpts: HashMap::with_capacity((est_fes * 3)/2 as uint), // estimate assumes most fes have 3 sides, each in 2 fes
      tags_by_fenum: if load_tags { Some(vec::with_capacity(est_fes)) } else { None },
      num_nb_sides: 0,
      num_b_sides: 0,
      sfs_btw_verts_work_set: HashSet::with_capacity(20),
      sf_endpt_pairs_work_buf: Some(vec::with_capacity(6)),
    };

    bldr.process_base_tris(base_tris, global_subdiv_iters);

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
                 bldr.take_tags())
  }


  fn process_base_tris
     ( &mut self,
       base_tris: &[BaseTri],
       global_subdiv_iters: uint )
  {
    let side_extra_subdiv_iters = side_extra_subdiv_iters_by_endpts(base_tris);

    // Process input mesh elements.
    for base_tri in base_tris.iter() {
      let nums_sfs = nums_side_faces_between_vertexes(base_tri, &side_extra_subdiv_iters);
      self.process_base_tri(base_tri, global_subdiv_iters, nums_sfs);
    }
  }
  

  fn process_base_tri
     ( &mut self,
       base_tri: &BaseTri,
       global_subdiv_iters: uint,
       nums_sfs_btw_verts: (u8, u8, u8) )
  {
    let (v0, v1, v2) = (base_tri.v0, base_tri.v1, base_tri.v2);
    
    // Add any extra subdivision iterations to be done in this element.
    let subdiv_iters = global_subdiv_iters + base_tri.extra_subdiv_iters();

    // Register the primary reference triangles for these subdivisions. Normally there will only be one such primary
    // (ie. non-inverted) reference triangle, which is just a scaled-down image of the base triangle. However when
    // multiple side faces are present between two vertexes, then additional reference triangles are required.
    let first_new_pri_oshape = self.register_primary_ref_tris(v0, v1, v2, nums_sfs_btw_verts, subdiv_iters);
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
                             base_tri.tag);
    } else { // no subdivision to be done
      // The base triangle itself is our finite element, with its own reference triangle.
      self.add_fe(first_new_pri_oshape, v0,v1,v2, base_tri.tag);
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
  fn register_primary_ref_tris
     ( &mut self,
       v0: Point, v1: Point, v2: Point,
       nums_sfs_btw_verts: (u8,u8,u8),
       subdiv_iters: uint )
     -> OShape
  {
    let first_new_oshape = OShape(self.oshapes.len());
   
    let scale = { let p: uint = pow_with_uint(2u, subdiv_iters); 1./(p as R) };

    if subdiv_iters == 0 || nums_sfs_btw_verts == (1u8, 1u8, 1u8) {
      self.oshapes.push(RefTri::new(v0,v1,v2, scale, nums_sfs_btw_verts));
    } else { 
      // We are subdividing, and there is more than one side face between some vertex pair. 

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
  fn register_secondary_ref_tri
     ( &mut self,
       pv0: Point, pv1: Point, pv2: Point, // base (primary) triangle vertexes
       subdiv_iters: uint ) 
     -> OShape
  {
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


  fn subdivide_primary
     ( &mut self,
       v0: Point, v1: Point, v2: Point,
       iters: uint,
       nums_side_faces_btw_verts: (u8,u8,u8),
       pri_oshapes_by_nums_sfs_btw_verts: |(u8,u8,u8)| -> OShape,
       sec_oshape: OShape,
       tag: Tag )
  {
    if iters == 0 {
      let oshape = pri_oshapes_by_nums_sfs_btw_verts(nums_side_faces_btw_verts);
      self.add_fe(oshape, v0,v1,v2, tag);
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
                             tag);

      self.subdivide_primary(midpt_v01, v1, midpt_v12,
                             iters-1,
                             (sfs_v01, sfs_v12, 1),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag);

      self.subdivide_primary(midpt_v20, midpt_v12, v2,
                             iters-1,
                             (1, sfs_v12, sfs_v20),
                             |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                             tag);

      // secondary sub-triangle
      self.subdivide_secondary(midpt_v01, midpt_v12, midpt_v20,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag);
    }
  }


  fn subdivide_secondary
     ( &mut self,
       v0: Point, v1: Point, v2: Point, // secondary triangle vertexes
       iters: uint,
       pri_oshapes_by_nums_sfs_btw_verts: |(u8,u8,u8)| -> OShape,
       sec_oshape: OShape,
       tag: Tag )
  {
    if iters == 0 { 
      self.add_fe(sec_oshape, v0,v1,v2, tag);
    } else {
      let (midpt_v01, midpt_v12, midpt_v20) = (midpt(v0,v1), midpt(v1,v2), midpt(v2,v0));

      // Corner subtriangles of this secondary triangle are secondary triangles needing one less subdivision iteration.
      self.subdivide_secondary(v0, midpt_v01, midpt_v20,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag);

      self.subdivide_secondary(midpt_v01, v1, midpt_v12,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag);

      self.subdivide_secondary(midpt_v20, midpt_v12, v2,
                               iters-1,
                               |sfs| pri_oshapes_by_nums_sfs_btw_verts(sfs), sec_oshape,
                               tag);

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
                             tag);
    }
  }


  // This procedure will be used to create all finite elements during the processing of the input base triangles.
  fn add_fe
     ( &mut self,
       oshape: OShape,
       v0: Point, v1: Point, v2: Point,
       tag: Tag )
  {
    // Create the finite element.
    let fe = FENum(self.fes.len());
    self.fes.push(ElTri{oshape: oshape, v0: v0});
    
    // optionally store tags
    match self.tags_by_fenum { Some(ref mut tags) => tags.push(tag), _ => {} }

    // Register this fe's local side representations ( fe/sf pairs ) by their endpoint pairs.
    let (nb_sides_delta, b_sides_delta) = {
      let nums_side_faces_btw_verts = self.oshapes[*oshape].nums_side_faces_between_vertexes;
      self.register_fe_side_reps_by_endpoints(fe, v0,v1,v2, nums_side_faces_btw_verts)
    };

    self.num_nb_sides = ((self.num_nb_sides as int) + nb_sides_delta) as uint;
    self.num_b_sides =  ((self.num_b_sides as int)  + b_sides_delta) as uint;
  }
 

  // Register the given finite element's side faces by their endpoints in the passed registry.
  // Returns the pair of changes in the number of non-boundary sides and boundary sides.
  fn register_fe_side_reps_by_endpoints
     ( &mut self, 
       fe: FENum,
       v0: Point, v1: Point, v2: Point,
       nums_side_faces_btw_verts: (u8,u8,u8) )
     -> (int, int)
  {
    let mut nb_sides_delta = 0;
    let mut b_sides_delta = 0;
    
    // Take the endpoint pairs work buffer, filled with the side face endpoint pairs.
    let endpt_pairs = TriMeshBuilder::fill_side_face_endpoint_pairs(self.take_side_face_endpoint_pairs_buf(),
                                                                    v0,v1,v2, 
                                                                    nums_side_faces_btw_verts,
                                                                    LesserEndpointsFirst);

    // Register the side face representations for this element by their endpoints.
    for sf in range(0, endpt_pairs.len()) {
      let fe_sf = (fe, SideFace(sf));
      self.side_reps_by_endpts.mangle(endpt_pairs[sf], (/*no context value*/),
        |_k, _ctx| { // no existing side rep found: insert
          b_sides_delta += 1;
          SideReps::new(fe_sf) // new value to insert
        },
        |_k, side_reps, _ctx| { // existing side rep found: mutate
          assert!(side_reps.snd_rep.is_none());
          side_reps.add(fe_sf);
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
  fn fill_side_face_endpoint_pairs
     ( mut endpts_buf: ~[(Point,Point)],
       v0: Point, v1: Point, v2: Point,
       nums_side_faces_btw_verts: (u8,u8,u8),
       endpoints_ordering: SideEndpointsOrdering )
     -> ~[(Point,Point)]
  {
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
  fn create_nb_sides_data
     ( &mut self )
     -> (StorageByInts2<Option<NBSideNum>>, ~[NBSideInclusions])
  {
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


  fn take_fes(&mut self) -> ~[ElTri] {
    let mut fes = ~[];
    swap(&mut fes, &mut self.fes);
    fes
  }
  
  fn take_oshapes(&mut self) -> ~[RefTri] {
    let mut oshapes = ~[];
    swap(&mut oshapes, &mut self.oshapes);
    oshapes
  }
  
  fn take_tags(&mut self) -> Option<~[Tag]> {
    let mut tags = None;
    swap(&mut tags, &mut self.tags_by_fenum);
    tags 
  }
  
  fn take_side_face_endpoint_pairs_buf(&mut self) -> ~[(Point,Point)] {
    let mut endpt_pairs_opt = None;
    swap(&mut endpt_pairs_opt, &mut self.sf_endpt_pairs_work_buf);
    endpt_pairs_opt.unwrap()
  }
  
  fn return_side_face_endpoint_pairs_buf(&mut self, buf: ~[(Point,Point)]) {
    let mut ret_buf_opt = Some(buf);
    swap(&mut ret_buf_opt, &mut self.sf_endpt_pairs_work_buf);
  }

} // TriMeshBuilder impl


// Represents an element in the input which may be subdivided to contribute elements to the final mesh.
struct BaseTri {

  v0: Point,
  v1: Point,
  v2: Point,
  
  tag: Tag,
}

impl BaseTri {
  
  fn extra_subdiv_iters(&self) -> uint {
    let Tag(subdivs) = self.tag;
    subdivs as uint
  }

} // BaseTri impl


// A structure to hold the one or two local representations of any side as fe/side face pairs.
struct SideReps {
  fst_rep: (FENum, SideFace),
  snd_rep: Option<(FENum, SideFace)>,
}

impl SideReps {
  
  fn new(fe_sf: (FENum,SideFace)) -> SideReps { SideReps { fst_rep: fe_sf, snd_rep: None } } 

  fn add(&mut self, (fe,sf): (FENum,SideFace)) {
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


// standalone utility functions


/* For each of the three sides of this base triangle to be subdivided, we determine here the number of faces
 * that the generated subtriangles should have between their vertexes which lie on the base triangle sides.
 * This allows these subdivision elements to meet those of a finer subdivision in a mesh element adjacent to this
 * base element ("hanging nodes"). The numbers are returned as a triplet of integers, corresponding to the side
 * face counts to be generated for elements' sides lying within sides v0v1, v1v2, and v2v0.
 */
fn nums_side_faces_between_vertexes <PTS_UINT_MAP: Map<(Point,Point),uint>>
   ( base_tri: &BaseTri,
     side_extra_subdiv_iters_by_endpts: &PTS_UINT_MAP )
   -> (u8, u8, u8)
{
  let base_tri_subdivs = base_tri.extra_subdiv_iters();
 
  let sfs_btw_verts = |va: Point, vb: Point| -> u8 {
    let sorted_endpts = mk_endpoint_pair(va, vb, LesserEndpointsFirst);
    let side_subdivs = match side_extra_subdiv_iters_by_endpts.find(&sorted_endpts) {
      Some(&subdivs) => subdivs,
      None => base_tri_subdivs
    };
    assert!(side_subdivs >= base_tri_subdivs);
    let sfs = pow_with_uint(2, side_subdivs - base_tri_subdivs);
    from_uint(sfs).unwrap()
  };

  (sfs_btw_verts(base_tri.v0, base_tri.v1),
   sfs_btw_verts(base_tri.v1, base_tri.v2),
   sfs_btw_verts(base_tri.v2, base_tri.v0))
}


/* 
 * Return a mapping of side endpoint pairs (lesser point first) to the the number of subdivision operations to be done
 * on the side in addition to the global number of subdivision operations. For a given non-boundary side this is the
 * larger of the numbers of extra subdivision operations for the two finite elements which include the side. No entry
 * will be present for sides for which the number of extra subdivisions is the same for both including finite elements.
 * Boundary sides may be included.
 */
fn side_extra_subdiv_iters_by_endpts
   ( base_tris: &[BaseTri] )
   -> HashMap<(Point,Point), uint>
{
  let mut side_subdivs_by_endpts: HashMap<(Point,Point),uint> = HashMap::new();

  for base_tri in base_tris.iter() {
    let base_tri_subdivs = base_tri.extra_subdiv_iters();
    let vert_pairs = [(base_tri.v0, base_tri.v1), (base_tri.v1, base_tri.v2), (base_tri.v2, base_tri.v0)];
    for &(va, vb) in vert_pairs.iter() {
      let sorted_endpts = mk_endpoint_pair(va, vb, LesserEndpointsFirst);
      let side_subdivs = side_subdivs_by_endpts.find_copy(&sorted_endpts);
      match side_subdivs {
        Some(subdivs) if subdivs == base_tri_subdivs => { // same # subdivs as prev, no entry needed
          side_subdivs_by_endpts.remove(&sorted_endpts);
        }
        Some(subdivs) if subdivs > base_tri_subdivs => {} // prev including base tri has more subdivs, leave as is
        _ => { // this base tri has more subdivs than prev including one, or this is the first occurence: insert
          side_subdivs_by_endpts.insert(sorted_endpts, base_tri_subdivs);
        }
      }
    }
  }

  side_subdivs_by_endpts
}


fn read_msh_nodes <I: Reader>
   ( msh_is: &mut BufferedReader<I> )
   -> ~[Point]
{
  // Position to the first line under the Nodes section.
  let mut lines = msh_is.lines().skip_while(|line| !line.starts_with("$Nodes")).skip(1);
  
  // First line in Nodes section should be the node count.
  let point_count: uint = match lines.next().map(|line| from_str::<uint>(line.trim_right())) {
    Some(Some(n)) => n,
    _ => fail!("Expected Nodes section beginning with node count - one of these was not found.")
  };
  
  let mut nodes = vec::from_elem(point_count, (R_NaN, R_NaN));  
 
  for line in lines.take_while(|line| !line.starts_with("$EndNodes")) {
    let node: Point = (extract_field(line, TOKIX_NODELINE_COORD1), extract_field(line, TOKIX_NODELINE_COORD2));
    let coord_3: Option<R> = maybe_extract_field(line, TOKIX_NODELINE_COORD3);
    match coord_3 {
      Some(nz) if nz != 0. => {
        fail!("Nodes with non-zero third coordinates are not supported in this 2D mesh reader.");
      }
      _ => {} // OK, no non-zero 3rd coordinate
    }
    let nodenum: uint = extract_field(line, TOKIX_NODELINE_NODENUM);
    nodes[nodenum-1] = node;
  }

  nodes
}


fn base_tri_from_line
   ( el_line:          &str,
     points_by_nodeix: &[Point] )
   -> Option<BaseTri>
{
  let el_type: uint = extract_field(el_line, TOKIX_ELLINE_ELTYPE);

  match el_type {
    ELTYPE_3_NODE_TRIANGLE => {
      let num_tags: uint = extract_field(el_line, TOKIX_ELLINE_NUMTAGS);
      assert!(num_tags >= 2);
      let last_tag_tokix = TOKIX_ELLINE_NUMTAGS + num_tags;

      let tag = Tag(extract_field(el_line, TOKIX_ELLINE_NUMTAGS+1));

      let nodenums = (extract_field::<uint>(el_line, last_tag_tokix + 1),
                      extract_field::<uint>(el_line, last_tag_tokix + 2),
                      extract_field::<uint>(el_line, last_tag_tokix + 3));
      
      let (v0, v1, v2) = (points_by_nodeix[nodenums.n0()-1],
                          points_by_nodeix[nodenums.n1()-1],
                          points_by_nodeix[nodenums.n2()-1]);
      
      Some(BaseTri{v0:v0, v1:v1, v2:v2, tag: tag})
    }
    el_type if is_lower_order_el_type(el_type) => None,
    _ => fail!("Element type not supported.")
  }
}


// Return the index of occurrence n of char c in string s, if any.
fn str_char_occur
   ( s: &str, c: char, n: uint )
   -> Option<uint>
{
  assert!(n > 0);
  let mut count = 0;
  for i in range(0, s.len()) {
    if s.char_at(i) == c {
      count += 1;
      if count == n { return Some(i); }
    }
  }
  None
}


#[inline]
fn extract_field <A: FromStr>
   ( line: &str, i: uint )
   -> A
{
  maybe_extract_field(line, i).unwrap()
}


#[inline]
fn maybe_extract_field <A: FromStr>
   ( line: &str, i: uint )
   -> Option<A>
{
  let start = if i == 0 { 0 } else { 1 + str_char_occur(line, ' ', i).unwrap() };
  let end = match str_char_occur(line, ' ', i+1) { Some(e) => e, None => line.len()-1 };
  from_str(line.slice(start, end))
}


fn is_lower_order_el_type
   ( el_type: uint )
   -> bool
{
  el_type == 15 || // point
  el_type == 1  || // 2 node line
  el_type == 8  || // 3 node line
  el_type == 26 || // 4 node line
  el_type == 27 || // 5 node line
  el_type == 28    // 6 node line
}


static ELTYPE_3_NODE_TRIANGLE: uint = 2;
static TOKIX_NODELINE_NODENUM: uint = 0;
static TOKIX_NODELINE_COORD1: uint = 1;
static TOKIX_NODELINE_COORD2: uint = 2;
static TOKIX_NODELINE_COORD3: uint = 3;
// Gmsh element line format: elm-number elm-type number-of-tags < tag > ... vert-number-list
static TOKIX_ELLINE_ELTYPE: uint = 1;
static TOKIX_ELLINE_NUMTAGS: uint = 2;


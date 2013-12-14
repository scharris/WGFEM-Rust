use common::{R, R_NaN};
use mesh::{Mesh, FENum, OShape, SideFace, NBSideNum, NBSideInclusions};
use wg_basis::{WGBasis, FaceMonNum};
use monomial::Monomial;
use tensor::{Tensor3, Tensor4, Tensor5};
use sparse_matrix::{SparseMatrix};

use std::vec;
use std::option::{Option};
use extra::sort;


pub trait VariationalBilinearForm<Mon:Monomial,MeshT:Mesh<Mon>> {

  fn basis<'a>(&'a self) -> &'a WGBasis<Mon,MeshT>;

  fn is_symmetric(&self) -> bool;

  fn int_mon_vs_int_mon(&self,
                        oshape: OShape,
                        monn_1: FaceMonNum,
                        monn_2: FaceMonNum) -> R;

  fn side_mon_vs_int_mon(&self,
                         oshape: OShape,
                         side_monn: FaceMonNum, side_face: SideFace,
                         int_monn: FaceMonNum) -> R;
  
  fn int_mon_vs_side_mon(&self,
                         oshape: OShape,
                         int_monn: FaceMonNum,
                         side_monn: FaceMonNum, side_face: SideFace) -> R;

  fn side_mon_vs_side_mon_fe_contr(&self,
                                   oshape: OShape,
                                   monn_1: FaceMonNum, side_face_1: SideFace,
                                   monn_2: FaceMonNum, side_face_2: SideFace) -> R;

  fn basis_els_vs_basis_els_transpose(&self) -> SparseMatrix {

    let (basis, mesh) = (self.basis(), self.basis().mesh());

    let sym = self.is_symmetric();
    let (num_int_mons, num_side_mons) = (basis.mons_per_fe_int(), basis.mons_per_fe_side());

    // The system matrix to be built.
    let mut m = SparseMatrix::new_with_capacities(basis.est_num_el_el_pairs_with_common_supp_fes(sym), basis.num_els());

    // Buffer to store non-boundary sides' numbers and fe side faces for either one or two finite elements.
    let mut fe_nb_sides_buf = vec::from_elem(mesh.max_num_shape_sides()*2, (NBSideNum(0), FENum(0), SideFace(0)));
    // Buffer to store the non-boundary side interactions with a single given non-boundary side.
    let mut nb_side_interactions_buf = vec::from_elem(mesh.max_num_shape_sides()*2,
      OneFESideSideInter(NBSideNum(0), NBSideNum(0), FENum(0), SideFace(0), SideFace(0)));

    // Create precomputed reference vbf values between monomials on single oriented shapes.
    let (int_vs_int_vbf_vals, int_vs_side_vbf_vals, side_vs_int_vbf_vals, side_vs_side_vbf_fe_contrs) = 
      (self.ref_int_vs_int_vbf_values(), 
       if !sym { self.ref_int_vs_side_vbf_values() } else { Tensor4::from_elem(0,0,0,0,0 as R) },
       self.ref_side_vs_int_vbf_values(),
       self.ref_side_vs_side_vbf_fe_contrs());

    // Iteration Order
    // We have to fill the sparse matrix in row-major order to satisfy the sparse matrix format requirments.
    // It's easiest to do this by iterating basis monomial pairs as they exist in the non-transposed matrix,
    // and then to insert the transposed values at the resulting row and column positions. In the case that
    // the vbf is symmetric, we only include pairs for the upper triangular part of the resulting matrix,
    // so that the second element number is greater than or equal to the first.
    // To iterate basis elements in order, it is necessary to iterate interior supported element first and
    // then side supported elements.  Within the interior supported elements, iteration procedes in order of
    // finite element number and interior monomial number, with monomial number less significant (changing
    // faster). For side supported elements, iteration is by increasing non-boundary side number and 
    // monomial number, with monomial number less significant.

    // Iterate basis element pairs beginning with an interior supported element. 
    
    for fe in range(0, mesh.num_fes()) { let fe = FENum(fe);
      let oshape = mesh.oriented_shape_for_fe(fe);

      for monn_1 in range(0, num_int_mons) { let monn_1 = FaceMonNum(monn_1);
        let r = *basis.int_mon_el_num(fe, monn_1);

        // Iterate element pairs with our interior supported element vs interior supported elements.
        for monn_2 in range(if sym { *monn_1 } else { 0 }, num_int_mons) { let monn_2 = FaceMonNum(monn_2);
          let c = *basis.int_mon_el_num(fe, monn_2);
          // As for all value fetches in this method, the value is fetched here with the second basis element in
          // first position because this is the *transpose* of the el vs el matrix.
          let ip = int_vs_int_vbf_vals.get(*oshape, *monn_2, *monn_1);
          if ip != 0. as R || r == c {
            m.push(r, c, ip);
          }
        }

        // Iterate element pairs with our interior supported element vs side supported elements.
        // The support sides are iterated in order of increasing non-boundary side number.
        for &(nbs, _, sf) in asc_nbsn_fe_sf_triplets_for_fes(fe, None, |_nbs_2| true, mesh, fe_nb_sides_buf).iter() {
          for monn_2 in range(0, num_side_mons) { let monn_2 = FaceMonNum(monn_2);
            let c = *basis.nb_side_mon_el_num(nbs, monn_2);
            let ip = side_vs_int_vbf_vals.get(*oshape, *monn_2, *sf, *monn_1);
            if ip != 0. as R {
              m.push(r, c, ip);
            }
          }
        }
      }
    }

    // Iterate basis element pairs beginning with a side supported element.

    for nbs in range(0, mesh.num_nb_sides()) { let nbs = NBSideNum(nbs);
      // Get the representations of our nbs as the side faces of finite elements.
      let nbs_incls = mesh.fe_inclusions_of_nb_side(nbs);

      let nbs_filter = if sym { |nbs_2: NBSideNum| nbs_2 >= nbs } else { |_nbs_2| true };
      let nb_side_interactions = asc_nb_sides_interactions(&nbs_incls,
                                                           nbs_filter,
                                                           mesh,
                                                           fe_nb_sides_buf, nb_side_interactions_buf);
  
      for monn_1 in range(0, num_side_mons) { let monn_1 = FaceMonNum(monn_1);
        let r = *basis.nb_side_mon_el_num(nbs, monn_1);
       
        // Iterate element pairs with our side supported element vs interior supported elements,
        // if the vbf is not symmetric.  A side vs interior element pairing is always in the
        // lower triangular part of the matrix which we won't provide in the symmetric case.
        if !sym {
          let incl_fe_sf_pairs = [(nbs_incls.fe1, nbs_incls.side_face_in_fe1), (nbs_incls.fe2, nbs_incls.side_face_in_fe2)];
          for &(mon_2_fe, mon_1_sf) in incl_fe_sf_pairs.iter() {
            let mon_2_fe_oshape = mesh.oriented_shape_for_fe(mon_2_fe);
            for monn_2 in range(0, num_int_mons) { let monn_2 = FaceMonNum(monn_2);
              let c = *basis.int_mon_el_num(mon_2_fe, monn_2);
              let ip = int_vs_side_vbf_vals.get(*mon_2_fe_oshape, *monn_2, *monn_1, *mon_1_sf);
              if ip != 0. as R {
                m.push(r, c, ip);
              }
            }
          }
        }
        
        // Iterate element pairs with our side supported element vs side supported elements.
       
        // Loop over the non-boundary side interactions with our first side, nbs.
        for nb_sides_interaction in nb_side_interactions.iter() {
          match nb_sides_interaction {
            &OneFESideSideInter(_, nbs_2, fe, nbs_sf_in_fe, nbs_2_sf_in_fe) => {
              let fe_oshape = mesh.oriented_shape_for_fe(fe);
              let monn_range_lower = if sym && nbs_2 == nbs { *monn_1 } else { 0 };
              for monn_2 in range(monn_range_lower, num_side_mons) { let monn_2 = FaceMonNum(monn_2);
                 let c = *basis.nb_side_mon_el_num(nbs_2, monn_2);
                 let ip = self.get_side_vs_side_vbf_contr(fe_oshape, monn_2, nbs_2_sf_in_fe, monn_1, nbs_sf_in_fe,
                                                          &side_vs_side_vbf_fe_contrs);
                 if ip != 0. as R || r == c {
                   m.push(r, c, ip);
                 }
              }
            }
            &TwoFESideSideInter(_, nbs_2,
                                fe_a, nbs_sf_in_fe_a, nbs_2_sf_in_fe_a,
                                fe_b, nbs_sf_in_fe_b, nbs_2_sf_in_fe_b) => {
              let fe_a_oshape = mesh.oriented_shape_for_fe(fe_a);
              let fe_b_oshape = mesh.oriented_shape_for_fe(fe_b);
              let monn_range_lower = if sym && nbs_2 == nbs { *monn_1 } else { 0 };
              for monn_2 in range(monn_range_lower, num_side_mons) { let monn_2 = FaceMonNum(monn_2);
                let c = *basis.nb_side_mon_el_num(nbs_2, monn_2);
                let ip = self.get_side_vs_side_vbf_contr(fe_a_oshape, monn_2, nbs_2_sf_in_fe_a, monn_1, nbs_sf_in_fe_a,
                                                         &side_vs_side_vbf_fe_contrs) + 
                         self.get_side_vs_side_vbf_contr(fe_b_oshape, monn_2, nbs_2_sf_in_fe_b, monn_1, nbs_sf_in_fe_b,
                                                         &side_vs_side_vbf_fe_contrs);
                if ip != 0. as R || r == c {
                  m.push(r, c, ip);
                }
              }
            }
          }
        } // side-side interactions
      } // monn_1
    } // nbs

    m
  }

  /* Returns a tensor of interior monomial vs interior monomial vbf values.  Results are indexed by oshape,
   * first monomial number, and second monomial number. If this variational form is symmetric, then only
   * values for which the first monomial is greater or equal to the second are provided.
   */
  fn ref_int_vs_int_vbf_values(&self) -> Tensor3 { // indexed by oshape, monn, monn
    let basis = self.basis();
    let num_int_mons = basis.mons_per_fe_int();
    let sym = self.is_symmetric();
    Tensor3::from_fn(basis.mesh().num_oriented_element_shapes(), num_int_mons, num_int_mons, |os, monn_1, monn_2| {
      if !sym || monn_1 >= monn_2 {
        self.int_mon_vs_int_mon(OShape(os), FaceMonNum(monn_1), FaceMonNum(monn_2))
      }
      else { R_NaN }
    })
  }

  /* Returns a tensor of vbf values indexed by oshape, side monomial number, side face, and interior monomial number.
   */
  fn ref_side_vs_int_vbf_values(&self) -> Tensor4 { // indexed by oshape, side monn, sf, int monn
    let basis = self.basis();
    Tensor4::from_fn(basis.mesh().num_oriented_element_shapes(),
                     basis.mons_per_fe_side(),
                     basis.mesh().max_num_shape_sides(),
                     basis.mons_per_fe_int(),
                     |os, side_monn, sf, int_monn| {
      if sf < basis.mesh().num_side_faces_for_oshape(OShape(os)) {
        self.side_mon_vs_int_mon(OShape(os), FaceMonNum(side_monn), SideFace(sf), FaceMonNum(int_monn))
      }
      else { R_NaN }
    })
  }

  /* Returns a tensor of vbf values indexed by oshape, interior monomial number, side monomial number, and side face.
   */
  fn ref_int_vs_side_vbf_values(&self) -> Tensor4 { // indexed by oshape, int monn, side monn, sf
    let basis = self.basis();
    Tensor4::from_fn(basis.mesh().num_oriented_element_shapes(),
                     basis.mons_per_fe_int(),
                     basis.mons_per_fe_side(),
                     basis.mesh().max_num_shape_sides(),
                     |os, int_monn, side_monn, sf| {
      if sf < basis.mesh().num_side_faces_for_oshape(OShape(os)) {
        self.int_mon_vs_side_mon(OShape(os), FaceMonNum(int_monn), FaceMonNum(side_monn), SideFace(sf))
      }
      else { R_NaN }
    })
  }

  /* Returns a tensor of contributions to side vs side vbf values on finite elements, indexed by finite element oshape,
   * first and second side faces, and first and second monomial numbers.  If the vbf is symmetric, then values are only
   * defined where either the first side face is greater than the second, or the side faces are equal and the first
   * monomial number is greater or equal to the second.  Callers should use the auxillary get_side_vs_side_vbf_contr
   * function to obtain values from the returned data, which will find the proper value by reversing indexes as
   * necessary in the case that the vbf is symmetric.
   */
  fn ref_side_vs_side_vbf_fe_contrs(&self) -> Tensor5 { // indexed by oshape, side monn, sf, side monn, sf
    let basis = self.basis();
    let num_oshapes = basis.mesh().num_oriented_element_shapes();
    let (max_sides, num_side_mons) = (basis.mesh().max_num_shape_sides(), basis.mons_per_fe_side());
    let sym = self.is_symmetric();
    Tensor5::from_fn(num_oshapes, num_side_mons, max_sides, num_side_mons, max_sides, |os, monn_1, sf_1, monn_2, sf_2| {
      let num_sides = basis.mesh().num_side_faces_for_oshape(OShape(os));
      if sf_1 >= num_sides || sf_2 >= num_sides || // non-existent side
         sym && (sf_1 < sf_2 || sf_1 == sf_2 && monn_1 < monn_2) { // value should be obtained via other indexes by symmetry
        R_NaN
      }
      else {
        self.side_mon_vs_side_mon_fe_contr(OShape(os), FaceMonNum(monn_1), SideFace(sf_1), FaceMonNum(monn_2), SideFace(sf_2))
      }
    })
  }

  #[inline]
  fn get_side_vs_side_vbf_contr(&self,
                                oshape: OShape, 
                                monn_1: FaceMonNum, sf_1: SideFace,
                                monn_2: FaceMonNum, sf_2: SideFace,
                                side_vs_side_vbf_fe_contrs: &Tensor5) -> R {
    if !self.is_symmetric() {
      side_vs_side_vbf_fe_contrs.get(*oshape, *monn_1, *sf_1, *monn_2, *sf_2)
    }
    else if sf_1 < sf_2 {
      side_vs_side_vbf_fe_contrs.get(*oshape, *monn_2, *sf_2, *monn_1, *sf_1)
    }
    else if sf_1 == sf_2 {
      let (lesser_monn, greater_monn) = if monn_1 < monn_2 { (monn_1, monn_2) } else { (monn_2, monn_1) };
      side_vs_side_vbf_fe_contrs.get(*oshape, *greater_monn, *sf_1, *lesser_monn, *sf_1)
    }
    else {
      side_vs_side_vbf_fe_contrs.get(*oshape, *monn_1, *sf_1, *monn_2, *sf_2)
    }
  }

} // trait VariationalBilinearForm


// Return (non-boundary side number, finite element, side face) triplets for all non-boundary sides of both of the given
// finite elements, in ascending order.
#[inline]
fn asc_nbsn_fe_sf_triplets_for_fes<'a, Mon:Monomial, MeshT:Mesh<Mon>>
   (fe1: FENum,
    fe2: Option<FENum>,
    nbs_filter: |NBSideNum| -> bool,
    mesh: &MeshT, fe_nb_sides_buf: &'a mut [(NBSideNum, FENum, SideFace)])
   -> &'a [(NBSideNum, FENum, SideFace)] {

  let mut next_pos = 0;

  let fe1_singleton = [fe1];
  for &fe in fe1_singleton.iter().chain(fe2.iter()) {
    for sf in range(0, mesh.num_side_faces_for_oshape(mesh.oriented_shape_for_fe(fe))) {
      if !mesh.is_boundary_side(fe, SideFace(sf)) {
        let nbs = mesh.nb_side_num_for_fe_side(fe, SideFace(sf));
        if nbs_filter(nbs) {
          fe_nb_sides_buf[next_pos] = (nbs, fe, SideFace(sf));
          next_pos += 1;
        }
      }
    }
  }

  let nb_sides = fe_nb_sides_buf.mut_slice(0, next_pos);
  sort::quick_sort3(nb_sides);
  nb_sides.as_slice()
}

// Returns all non-boundary side interactions for the side with the given finite element inclusions, in order of increasing
// second non-boundary side numbers (the first being fixed).
#[inline]
fn asc_nb_sides_interactions<'a, Mon:Monomial, MeshT:Mesh<Mon>>
   (nbs_1_incls: &NBSideInclusions,
    nbs_filter: |NBSideNum| -> bool,
    mesh: &MeshT,
    fe_nb_sides_buf: &mut [(NBSideNum, FENum, SideFace)],
    nb_side_interactions_buf: &'a mut [NBSidesInteraction]) -> &'a [NBSidesInteraction] {

  let nbs_1 = nbs_1_incls.nb_side_num;
  let mut nbs_2_it = asc_nbsn_fe_sf_triplets_for_fes(nbs_1_incls.fe1, Some(nbs_1_incls.fe2), nbs_filter, mesh, fe_nb_sides_buf).iter().peekable();
  
  let mut next_inter = 0u;

  loop {
    match nbs_2_it.next() {
      None => { break; },
      Some(&(nbs_2, fe_a, nbs_2_sf_in_fe_a)) => {
        let nbs_1_sf_in_fe_a = match fe_a { 
          fe if fe == nbs_1_incls.fe1 => nbs_1_incls.side_face_in_fe1,
          fe if fe == nbs_1_incls.fe2 => nbs_1_incls.side_face_in_fe2,
          _ => fail!("Unexpected finite element: expected one of including fe's of first non-boundary side.")
        };
        match nbs_2_it.peek() {
          Some(& &(some_nbs, _, _)) if some_nbs == nbs_2 => { // Side nbs_2 is also included in the other of nbs_1's including fe's.
            let (fe_b, nbs_2_sf_in_fe_b) = match nbs_2_it.next() { Some(&(_, fe, sf)) => (fe, sf),
                                                                   None => fail!("Peeked value could not be obtained(?)") };
            let nbs_1_sf_in_fe_b = match fe_b { 
              fe if fe == nbs_1_incls.fe1 => nbs_1_incls.side_face_in_fe1,
              fe if fe == nbs_1_incls.fe2 => nbs_1_incls.side_face_in_fe2,
              _ => fail!("Unexpected finite element: expected one of including fe's of first non-boundary side.")
            };

            nb_side_interactions_buf[next_inter] = TwoFESideSideInter(nbs_1, nbs_2,
                                                                      fe_a, nbs_1_sf_in_fe_a, nbs_2_sf_in_fe_a,
                                                                      fe_b, nbs_1_sf_in_fe_b, nbs_2_sf_in_fe_b);
            next_inter += 1;
          }
          _ => {
            nb_side_interactions_buf[next_inter] = OneFESideSideInter(nbs_1, nbs_2, fe_a, nbs_1_sf_in_fe_a, nbs_2_sf_in_fe_a);
            next_inter += 1;
          }
        }
      }
    }
  }

  nb_side_interactions_buf.slice(0, next_inter)
}

/* NBSidesInteraction
 * Represents the interaction of a pair of non-boundary sides having at least one common support finite element.
 * The interaction is captured for each common supporting finite element as the pair of side faces that the 
 * non-boundary sides occupy in the finite element.
 *
 * Sides which interact on a single finite element.
 *      ____
 *     |    | 
 *     |_s1_|_____
 *     |    |     |
 *     |    |s2   | 
 *  fe |____|_____|
 *
 *
 * Sides which interact on two finite elements, fe_a and fe_b.  In the first figure we could have s1 = s2 or also
 * s1 and s2 representing two sides with a straight angle between them ("hanging node").
 *         _________                ________
 *        |    |    |              |        | 
 *        |    |s1  |              |_s1_    |
 *        |    |    |              |    |   |
 *        |    |s2  |              |    |s2 |
 *   fe_a |____|____| fe_b    fe_a |____|___| fe_b
 */
#[deriving(Clone)]
enum NBSidesInteraction {

  OneFESideSideInter(NBSideNum, NBSideNum, FENum, SideFace, SideFace),

  TwoFESideSideInter(NBSideNum, NBSideNum, 
                     FENum, SideFace, SideFace,
                     FENum, SideFace, SideFace)
}


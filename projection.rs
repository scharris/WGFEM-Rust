use common::{R};
use la;
use la::lapack_int;
use dense_matrix::DenseMatrix;
use monomial::{Monomial};
use polynomial::{PolyBorrowingMons, PolyOwning};
use mesh::{Mesh, FENum, OShape, SideFace};
use wg_basis::{WGBasis};

use std::num;
use std::vec;

/** Method of Projection onto a Subspace
 * -------------------------------------
 * <Proj(g), e_i> = <g, e_i> for all basis elements e_i in the subspace.
 * Since Proj(g) = sum_j{ a_j e_j } for some a_j's, this gives us a linear system
 *   sum_j { <e_j, e_i> a_j } = <g, e_i>
 * which we solve for the a_j.
 */

pub struct Projector<'a,Mon,MeshT> {
 
  basis: &'a WGBasis<Mon,MeshT>,

  // Work matrices for lapack to avoid allocations and because lapack enjoys writing over its inputs.
  la_ips_int_mons: DenseMatrix,
  la_ips_side_mons: DenseMatrix,
  la_pivots: ~[lapack_int],
  la_pivots_buf: *mut lapack_int,
  la_int_proj_rhs: DenseMatrix,
  la_side_proj_rhs: DenseMatrix,
}

impl <'a,Mon:Monomial,MeshT:Mesh<Mon>> Projector<'a,Mon,MeshT> {

  pub fn new(basis: &'a WGBasis<Mon,MeshT>) -> Projector<'a,Mon,MeshT> {
    Projector::with_rhs_cols_capacity(basis, 1000u)
  }
  
  pub fn with_rhs_cols_capacity(basis: &'a WGBasis<Mon,MeshT>, init_rhs_cols_capacity: uint) -> Projector<'a,Mon,MeshT> {
    let (num_int_mons, num_side_mons) = (basis.mons_per_fe_int(), basis.mons_per_fe_side());
    let la_ips_int_mons = DenseMatrix::from_elem(num_int_mons, num_int_mons, 0 as R);  
    let la_ips_side_mons = DenseMatrix::from_elem(num_side_mons, num_side_mons, 0 as R);  
    let mut la_pivots = vec::from_elem(num::max(num_int_mons, num_side_mons), 0 as lapack_int);
    let la_pivots_buf = la_pivots.as_mut_ptr();
    let init_num_rhs_cols = init_rhs_cols_capacity;
    let la_int_proj_rhs = DenseMatrix::from_elem(num_int_mons, init_num_rhs_cols, 0 as R);
    let la_side_proj_rhs = DenseMatrix::from_elem(num_side_mons, init_num_rhs_cols, 0 as R);

    Projector {
      basis: basis,
      la_ips_int_mons: la_ips_int_mons,
      la_ips_side_mons: la_ips_side_mons,
      la_pivots: la_pivots,
      la_pivots_buf: la_pivots_buf,
      la_int_proj_rhs: la_int_proj_rhs,
      la_side_proj_rhs: la_side_proj_rhs,
    }
  }
 
  #[inline]
  pub fn basis(&self) -> &'a WGBasis<Mon,MeshT> {
    self.basis
  }

  
  /// Returns the projections of function g to the spaces spanned by the basis functions supported on
  /// the interiors of the given finite elements in turn. The monomials of the returned polynomials are
  /// gauranteed to be in order of their ascending face monomial numbers, which is the same order in
  /// which the monomials appear in the basis.
  #[inline(never)]
  pub fn projs_to_int_supp_approx_spaces(&mut self,
                                         g: |&[R]| -> R,
                                         fes: &[FENum],
                                         fes_oshape: OShape) -> ~[PolyBorrowingMons<'a,Mon>] {

    let int_mons = self.basis.ref_int_mons();

    let sol_as_col_maj_vec = unsafe {

      // Prepare system matrix a.
      let a = {
        self.la_ips_int_mons.fill_upper_triangle_from(self.basis.ips_int_mons_for_oshape(fes_oshape));
        self.la_ips_int_mons.mut_col_maj_data_ptr()
      };

      // Prepare right hand side b matrix, one column per fe to be projected onto.
      let b = {
        //   - Set the proper number of columns for our sequence of fes, recreate matrix if necessary.
        if self.la_int_proj_rhs.capacity_cols() >= fes.len() {
          self.la_int_proj_rhs.set_num_cols(fes.len());
        } else {
          self.la_int_proj_rhs = DenseMatrix::of_size_with_cols_capacity(int_mons.len(), fes.len(), 3*fes.len()/2);
        }
        //   - Fill the rhs matrix with inner products between g and our interior monomials (rows) on the fes (cols). 
        self.la_int_proj_rhs.fill_from_fn(|r,c| {
          let mon = int_mons[r].clone();
          let fe = fes[c]; 
          self.basis.mesh.intg_global_fn_x_facerel_mon_on_fe_int(|x| g(x), mon, fe)
        });
        self.la_int_proj_rhs.mut_col_maj_data_ptr()
      };

      la::solve_symmetric_as_col_maj_with_ut_sys(a, int_mons.len() as lapack_int,
                                                 b, fes.len() as lapack_int,
                                                 self.la_pivots_buf);

      vec::from_buf(b as *R, int_mons.len() * fes.len())
    };

    sol_as_col_maj_vec
      .chunks(int_mons.len()) // Chunk into columns representing the projection coefficients.
      .map(|proj_coefs| PolyBorrowingMons { coefs: proj_coefs.to_owned(), mons: int_mons })
      .collect()
  }
  
  /// Returns the projections of function g to the spaces spanned by the basis functions supported on
  /// the given side face of the given finite elements in turn. The monomials of the returned polynomials
  /// are gauranteed to be in order of their ascending face monomial numbers, which is the same order in
  /// which the monomials appear in the basis.
  #[inline(never)]
  pub fn projs_to_side_supp_approx_spaces(&mut self,
                                          g: |&[R]| -> R,
                                          fes: &[FENum],
                                          fes_oshape: OShape,
                                          side_face: SideFace) -> ~[PolyBorrowingMons<'a,Mon>] {

    let side_mons = self.basis.side_mons_for_oshape_side(fes_oshape, side_face);

    let sol_as_col_maj_vec = unsafe {

      // Prepare system matrix a.
      let a = {
        self.la_ips_side_mons.fill_upper_triangle_from(self.basis.ips_side_mons_for_oshape_side(fes_oshape, side_face));
        self.la_ips_side_mons.mut_col_maj_data_ptr()
      };

      // Prepare right hand side b matrix, one column per fe to be projected onto.
      let b = {
        //   - Set the proper number of columns for our sequence of fes, recreate matrix if necessary.
        if self.la_side_proj_rhs.capacity_cols() >= fes.len() {
          self.la_side_proj_rhs.set_num_cols(fes.len());
        } else {
          self.la_side_proj_rhs = DenseMatrix::of_size_with_cols_capacity(side_mons.len(), fes.len(), 3*fes.len()/2);
        }
        //   - Fill the rhs matrix with inner products between g and our side monomials (rows) on the fes (cols). 
        self.la_side_proj_rhs.fill_from_fn(|r,c| {
          let mon = side_mons[r].clone();
          let fe = fes[c]; 
          self.basis.mesh.intg_global_fn_x_facerel_mon_on_fe_side(|x| g(x), mon, fe, side_face)
        });
        self.la_side_proj_rhs.mut_col_maj_data_ptr()
      };

      la::solve_symmetric_as_col_maj_with_ut_sys(a, side_mons.len() as lapack_int,
                                                 b, fes.len() as lapack_int,
                                                 self.la_pivots_buf);

      vec::from_buf(b as *R, side_mons.len() * fes.len())
    };

    sol_as_col_maj_vec
      .chunks(side_mons.len()) // Chunk into columns representing the projection coefficients.
      .map(|proj_coefs| PolyBorrowingMons { coefs: proj_coefs.to_owned(), mons: side_mons })
      .collect()
  }

  #[inline(never)]
  pub fn projs_int_mons_to_side_supp_approx_space(&mut self,
                                                  int_mons: &[Mon],
                                                  oshape: OShape,
                                                  side_face: SideFace) -> ~[PolyOwning<Mon>] {

    let side_mons = self.basis.side_mons_for_oshape_side(oshape, side_face);
    
    let sol_as_col_maj_vec = unsafe {

      // Prepare system matrix a, consisting of inner products between basis elements supported on the side.
      let a = {
        self.la_ips_side_mons.fill_upper_triangle_from(self.basis.ips_side_mons_for_oshape_side(oshape, side_face));
        self.la_ips_side_mons.mut_col_maj_data_ptr()
      };

      // Prepare right hand side b matrix, one column per interior monomial to be projected onto the side.
      let b = {
        //   - Set the proper number of columns for our sequence of fes, recreate matrix if necessary.
        if self.la_side_proj_rhs.capacity_cols() >= int_mons.len() {
          self.la_side_proj_rhs.set_num_cols(int_mons.len());
        } else {
          self.la_side_proj_rhs = DenseMatrix::of_size_with_cols_capacity(side_mons.len(), int_mons.len(), 3*int_mons.len()/2);
        }
        //   - Fill the rhs matrix with inner products between our side monomials (rows) and the interior monomials (cols) to be projected. 
        self.la_side_proj_rhs.fill_from_fn(|r,c| {
          let (int_mon, side_basis_mon) = (int_mons[c].clone(), side_mons[r].clone());
          self.basis.mesh.intg_intrel_mon_x_siderel_mon_on_oshape_side(int_mon, side_basis_mon, oshape, side_face)
        });
        self.la_side_proj_rhs.mut_col_maj_data_ptr()
      };

      la::solve_symmetric_as_col_maj_with_ut_sys(a, side_mons.len() as lapack_int,
                                                 b, int_mons.len() as lapack_int,
                                                 self.la_pivots_buf);

      vec::from_buf(b as *R, side_mons.len() * int_mons.len())
    };

    sol_as_col_maj_vec
      .chunks(side_mons.len()) // Chunk into columns representing the projection coefficients.
      .map(|proj_coefs| PolyOwning { coefs: proj_coefs.to_owned(), mons: side_mons.to_owned() })
      .collect()
  }

} // impl Projector


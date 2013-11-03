use common::{R};
use lapack;
use lapack::lapack_int;
use dense_matrix::DenseMatrix;
use monomial::{Monomial};
use polynomial::{PolyBorrowingMons};
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

pub struct Projector<'self,Mon,MeshT> {
 
  basis: &'self WGBasis<Mon,MeshT>,

  // Work matrices for lapack to avoid allocations and because lapack enjoys writing over its inputs.
  lapack_ips_int_mons: DenseMatrix,
  lapack_ips_side_mons: DenseMatrix,
  lapack_pivots: ~[lapack_int],
  lapack_pivots_buf: *mut lapack_int,
  lapack_int_proj_rhs: DenseMatrix,
  lapack_side_proj_rhs: DenseMatrix,
}

impl <'self,Mon:Monomial,MeshT:Mesh<Mon>> Projector<'self,Mon,MeshT> {

  pub fn new(basis: &'self WGBasis<Mon,MeshT>) -> Projector<'self,Mon,MeshT> {
    Projector::with_rhs_cols_capacity(basis, 1000u)
  }
  
  pub fn with_rhs_cols_capacity(basis: &'self WGBasis<Mon,MeshT>, init_rhs_cols_capacity: uint) -> Projector<'self,Mon,MeshT> {
    lapack::init(); // TODO: Move this to main when available.

    let (num_int_mons, num_side_mons) = (basis.mons_per_fe_int(), basis.mons_per_fe_side());
    let lapack_ips_int_mons = DenseMatrix::from_elem(num_int_mons, num_int_mons, 0 as R);  
    let lapack_ips_side_mons = DenseMatrix::from_elem(num_side_mons, num_side_mons, 0 as R);  
    let mut lapack_pivots = vec::from_elem(num::max(num_int_mons, num_side_mons), 0 as lapack_int);
    let lapack_pivots_buf = vec::raw::to_mut_ptr(lapack_pivots);
    let init_num_rhs_cols = init_rhs_cols_capacity;
    let lapack_int_proj_rhs = DenseMatrix::from_elem(num_int_mons, init_num_rhs_cols, 0 as R);
    let lapack_side_proj_rhs = DenseMatrix::from_elem(num_side_mons, init_num_rhs_cols, 0 as R);

    Projector {
      basis: basis,
      lapack_ips_int_mons: lapack_ips_int_mons,
      lapack_ips_side_mons: lapack_ips_side_mons,
      lapack_pivots: lapack_pivots,
      lapack_pivots_buf: lapack_pivots_buf,
      lapack_int_proj_rhs: lapack_int_proj_rhs,
      lapack_side_proj_rhs: lapack_side_proj_rhs,
    }
  }
 
  #[inline]
  pub fn basis(&self) -> &'self WGBasis<Mon,MeshT> {
    self.basis
  }

  
  /// Returns the projections of function g to the spaces spanned by the basis functions supported on
  /// the interiors of the given finite elements in turn. The monomials of the returned polynomials are
  /// gauranteed to be in order of their ascending face monomial numbers, which is the same order in
  /// which the monomials appear in the basis.
  #[fixed_stack_segment]
  #[inline(never)]
  pub fn projs_to_span_fes_int_supp_basis_els(&mut self,
         g: &fn(&[R]) -> R, fes: &[FENum], fes_oshape: OShape) -> ~[PolyBorrowingMons<'self,Mon>] {

    let int_mons = self.basis.ref_int_mons();

    let sol_as_col_maj_vec = unsafe {

      // Prepare system matrix a.
      let a = {
        self.lapack_ips_int_mons.fill_upper_triangle_from(self.basis.ips_int_mons_for_oshape(fes_oshape));
        self.lapack_ips_int_mons.mut_col_maj_data_ptr()
      };

      // Prepare right hand side b matrix, one column per fe to be projected onto.
      let b = {
        //   - Set the proper number of columns for our sequence of fes, recreate matrix if necessary.
        if self.lapack_int_proj_rhs.capacity_cols >= fes.len() {
          self.lapack_int_proj_rhs.set_num_cols(fes.len());
        } else {
          self.lapack_int_proj_rhs = DenseMatrix::of_size_with_cols_capacity(int_mons.len(), fes.len(), 3*fes.len()/2);
        }
        //   - Fill the rhs matrix with inner products between g and our interior monomials (rows) on the fes (cols). 
        self.lapack_int_proj_rhs.fill_from_fn(|r,c| {
          let mon = int_mons[r].clone();
          let fe = fes[c]; 
          self.basis.mesh.intg_global_fn_x_facerel_mon_on_fe_int(|x| g(x), mon, fe)
        });
        self.lapack_int_proj_rhs.mut_col_maj_data_ptr()
      };

      lapack::solve_symmetric_as_col_maj_with_ut_sys(a, int_mons.len() as lapack_int,
                                                     b, fes.len() as lapack_int,
                                                     self.lapack_pivots_buf);

      vec::from_buf(b as *R, int_mons.len() * fes.len())
    };

    sol_as_col_maj_vec
      .chunk_iter(int_mons.len()) // Chunk into columns representing the projection coefficients.
      .map(|proj_coefs| PolyBorrowingMons { coefs: proj_coefs.to_owned(), mons: int_mons })
      .collect()
  }
  
  /// Returns the projections of function g to the spaces spanned by the basis functions supported on
  /// the given side face of the given finite elements in turn. The monomials of the returned polynomials
  /// are gauranteed to be in order of their ascending face monomial numbers, which is the same order in
  /// which the monomials appear in the basis.
  #[fixed_stack_segment]
  #[inline(never)]
  pub fn projs_to_span_fes_side_supp_basis_els(&mut self,
         g: &fn(&[R]) -> R, fes: &[FENum], fes_oshape: OShape, side_face: SideFace) -> ~[PolyBorrowingMons<'self,Mon>] {

    let side_mons = self.basis.side_mons_for_oshape_side(fes_oshape, side_face);

    let sol_as_col_maj_vec = unsafe {

      // Prepare system matrix a.
      let a = {
        self.lapack_ips_side_mons.fill_upper_triangle_from(self.basis.ips_side_mons_for_oshape_side(fes_oshape, side_face));
        self.lapack_ips_side_mons.mut_col_maj_data_ptr()
      };

      // Prepare right hand side b matrix, one column per fe to be projected onto.
      let b = {
        //   - Set the proper number of columns for our sequence of fes, recreate matrix if necessary.
        if self.lapack_side_proj_rhs.capacity_cols >= fes.len() {
          self.lapack_side_proj_rhs.set_num_cols(fes.len());
        } else {
          self.lapack_side_proj_rhs = DenseMatrix::of_size_with_cols_capacity(side_mons.len(), fes.len(), 3*fes.len()/2);
        }
        //   - Fill the rhs matrix with inner products between g and our side monomials (rows) on the fes (cols). 
        self.lapack_side_proj_rhs.fill_from_fn(|r,c| {
          let mon = side_mons[r].clone();
          let fe = fes[c]; 
          self.basis.mesh.intg_global_fn_x_facerel_mon_on_fe_side(|x| g(x), mon, fe, side_face)
        });
        self.lapack_side_proj_rhs.mut_col_maj_data_ptr()
      };

      lapack::solve_symmetric_as_col_maj_with_ut_sys(a, side_mons.len() as lapack_int,
                                                     b, fes.len() as lapack_int,
                                                     self.lapack_pivots_buf);

      vec::from_buf(b as *R, side_mons.len() * fes.len())
    };

    sol_as_col_maj_vec
      .chunk_iter(side_mons.len()) // Chunk into columns representing the projection coefficients.
      .map(|proj_coefs| PolyBorrowingMons { coefs: proj_coefs.to_owned(), mons: side_mons })
      .collect()
  }

  #[fixed_stack_segment]
  #[inline(never)]
  pub fn proj_int_mons_to_span_oshape_side_supp_basis_els(&mut self,
     int_mons: &[Mon], oshape: OShape, side_face: SideFace) -> ~[PolyBorrowingMons<'self,Mon>] {

    let side_mons = self.basis.side_mons_for_oshape_side(oshape, side_face);
    
    let sol_as_col_maj_vec = unsafe {

      // Prepare system matrix a, consisting of inner products between basis elements supported on the side.
      let a = {
        self.lapack_ips_side_mons.fill_upper_triangle_from(self.basis.ips_side_mons_for_oshape_side(oshape, side_face));
        self.lapack_ips_side_mons.mut_col_maj_data_ptr()
      };

      // Prepare right hand side b matrix, one column per interior monomial to be projected onto the side.
      let b = {
        //   - Set the proper number of columns for our sequence of fes, recreate matrix if necessary.
        if self.lapack_side_proj_rhs.capacity_cols >= int_mons.len() {
          self.lapack_side_proj_rhs.set_num_cols(int_mons.len());
        } else {
          self.lapack_side_proj_rhs = DenseMatrix::of_size_with_cols_capacity(side_mons.len(), int_mons.len(), 3*int_mons.len()/2);
        }
        //   - Fill the rhs matrix with inner products between our side monomials (rows) and the interior monomials (cols) to be projected. 
        self.lapack_side_proj_rhs.fill_from_fn(|r,c| {
          let (int_mon, side_basis_mon) = (int_mons[c].clone(), side_mons[r].clone());
          self.basis.mesh.intg_intrel_mon_x_siderel_mon_on_oshape_side(int_mon, side_basis_mon, oshape, side_face)
        });
        self.lapack_side_proj_rhs.mut_col_maj_data_ptr()
      };

      lapack::solve_symmetric_as_col_maj_with_ut_sys(a, side_mons.len() as lapack_int,
                                                     b, int_mons.len() as lapack_int,
                                                     self.lapack_pivots_buf);

      vec::from_buf(b as *R, side_mons.len() * int_mons.len())
    };

    sol_as_col_maj_vec
      .chunk_iter(side_mons.len()) // Chunk into columns representing the projection coefficients.
      .map(|proj_coefs| PolyBorrowingMons { coefs: proj_coefs.to_owned(), mons: side_mons })
      .collect()
  }

} // impl Projector


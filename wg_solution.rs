use common::R;
use monomial::Monomial;
use polynomial::{PolyBorrowing, PolyBorrowingMons};
use mesh::{Mesh, FENum, SideFace};
use wg_basis::{WGBasis, FaceMonNum};

use std::hashmap::HashMap;


pub type BoundaryProjections<'self,Mon> = HashMap<(FENum,SideFace), PolyBorrowingMons<'self,Mon>>;

pub struct WGSolution<'self,Mon,MeshT> {
  basis_coefs: ~[R],
  basis: &'self WGBasis<Mon,MeshT>,
  bnd_projs: BoundaryProjections<'self,Mon>
}

impl<'self,Mon:Monomial,MeshT:Mesh<Mon>> WGSolution<'self,Mon,MeshT> {

  pub fn new(basis_coefs: ~[R], basis: &'self WGBasis<Mon,MeshT>, bnd_projs: BoundaryProjections<'self,Mon>) -> WGSolution<'self,Mon,MeshT> {
    WGSolution {
      basis_coefs: basis_coefs,
      basis: basis,
      bnd_projs: bnd_projs,
    }
  }
  
  #[inline]
  pub fn basis(&self) -> &'self WGBasis<Mon,MeshT> {
    self.basis
  }
  
  #[inline]
  pub fn basis_coefs<'a>(&'a self) -> &'a [R] {
    self.basis_coefs.as_slice()
  }

  #[inline]
  pub fn bnd_projs<'a>(&'a self) -> &'a BoundaryProjections<'self,Mon> {
    &self.bnd_projs
  }

  #[inline]
  pub fn value_at_int_rel(&self, fe: FENum, x: &[R]) -> R {
    let fe_first_int_beln = self.basis.int_mon_el_num(fe, FaceMonNum(0));
    self.basis_coefs.slice_from(*fe_first_int_beln).iter()
                    .zip(self.basis.ref_int_mons().iter())
                    .fold(0 as R, |sum, (&coef, mon)| sum + coef * mon.value_at(x))
  }

  /// Get the polynomial representing the passed full WG solution restricted to a particular finite element interior.
  #[inline]
  pub fn fe_int_poly<'a>(&'a self, fe: FENum) -> PolyBorrowing<'a,Mon> {
    self.basis.fe_int_poly(fe, self.basis_coefs)
  }

  /// Get the polynomial representing the passed full WG solution restricted to a particular finite element interior.
  #[inline]
  pub fn fe_side_poly<'a>(&'a self, fe: FENum, side_face: SideFace) -> PolyBorrowing<'a,Mon> {
    self.basis.fe_side_poly(fe, side_face, self.basis_coefs)
  }

}



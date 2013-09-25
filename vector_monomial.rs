use std::vec;
use common::*;
use monomial::{Monomial};

 /*
  * VectorMonomial type, parameterized by monomial type M.
  * Represents a vector valued function with a monomial in one component and
  * which is constantly 0 in all other components.
  */
#[deriving(Eq, IterBytes, Clone)]
pub struct VectorMonomial<M> {
  mon_dim: Dim,
  mon: M,
}

impl<M:Monomial> VectorMonomial<M> {

  pub fn new(mon_dim: Dim, mon: M) -> VectorMonomial<M> {
    if *mon_dim >= Monomial::domain_space_dims(None::<M>) {
      fail!("Coordinate index too large for passed monomial's domain.");
    }
    VectorMonomial { mon_dim: mon_dim, mon: mon }
  }

  #[inline]
  pub fn mon_dim(&self) -> Dim {
    self.mon_dim
  }
  
  #[inline]
  pub fn mon(&self) -> M {
    self.mon.clone()
  }

  #[inline]
  pub fn divergence_coef_and_mon(&self) -> (R,M) {
    match self.mon.exp(self.mon_dim) {
      Deg(0) => (0.,self.mon.clone()),
      Deg(e) => (e as R, self.mon.map_exp(self.mon_dim, |e| Deg(*e-1)))
    }
  }

  /// Vector monomials of degree not exceeding that indicated, ordered by component dimension in ascending order,
  /// and then by monomial in increasing lexicographical order of exponents.
  pub fn vector_mons_of_deg_le(max_deg: Deg) -> ~[VectorMonomial<M>] {
    let dom_dim = Monomial::domain_space_dims(None::<M>);
    let mons: ~[M] = Monomial::mons_of_deg_le(max_deg); 
    let mut vmons: ~[VectorMonomial<M>] = vec::with_capacity(mons.len() * dom_dim);
    for r in range(0, dom_dim) {
      for mon in mons.iter() {
        vmons.push(VectorMonomial::new(Dim(r), mon.clone()));
      }
    }
    vmons
  }
}


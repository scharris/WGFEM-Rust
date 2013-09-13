use common::*;
use monomial::*;

mod common;
mod monomial;

// VectorMonomial, parameterized by monomial type M.
// Represents a vector valued function with a monomial in one component and
// which is constantly 0 in all other components.
#[deriving(Eq,IterBytes,Clone)]
pub struct VectorMonomial<M> {
  mon_coord_ix: Dim,
  mon: M,
}

impl<M:Monomial> VectorMonomial<M> {

  fn new(mon_coord_ix: Dim, mon: M) -> VectorMonomial<M> {
    if *mon_coord_ix >= *Monomial::domain_dim(None::<M>) {
      fail!("Coordinate index too large for passed monomial's domain.");
    }
    VectorMonomial { mon_coord_ix: mon_coord_ix, mon: mon }
  }

}


#[test]
fn test_construction() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x);
}

#[test]
#[should_fail]
fn test_improper_construction() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(3), x);
}



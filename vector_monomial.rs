use common::*;
use monomial::*;

// VectorMonomial, parameterized by monomial type M.
// Represents a vector valued function with a monomial in one component and
// which is constantly 0 in all other components.
#[deriving(Eq, IterBytes, Clone)]
pub struct VectorMonomial<M> {
  mon_dim: Dim,
  mon: M,
}

impl<M:Monomial> VectorMonomial<M> {

  pub fn new(mon_dim: Dim, mon: M) -> VectorMonomial<M> {
    if *mon_dim >= Monomial::domain_dim(None::<M>) {
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

#[test]
fn test_divergence() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  
  let x_dim0_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x);
  assert_eq!(x_dim0_vmon.divergence_coef_and_mon(), (1.,one));
  
  let y_dim0_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), y);
  assert_eq!(y_dim0_vmon.divergence_coef_and_mon(), (0.,y));
  
  let y_dim1_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), y);
  assert_eq!(y_dim1_vmon.divergence_coef_and_mon(), (1.,one));
  
  let xy_dim0_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x*y);
  assert_eq!(xy_dim0_vmon.divergence_coef_and_mon(), (1.,y));

  let x2y3_dim1_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x*x*y*y*y);
  assert_eq!(x2y3_dim1_vmon.divergence_coef_and_mon(), (3., x*x*y*y));
}


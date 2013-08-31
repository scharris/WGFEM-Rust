extern mod extra;
use extra::treemap::TreeMap;
use std::vec;
use common::*;
use monomials::*;

mod common;
mod monomials;

// Polynomial, parameterized by monomial type M
#[deriving(Clone)]
pub struct Polynomial<M> {
  // TODO: Maybe store as ~[(R,M)] instead
  mons: ~[M],
  coefs: ~[R]
}

impl <M:Monomial> Polynomial<M> {

  fn new(mons: ~[M], coefs: ~[R]) -> Polynomial<M> {
    match mons.len() {
      len if len != coefs.len() => { fail!("size of coefficent and monomial arrays must match"); }
      0 =>  { fail!("polynomial requires one or more terms"); }
      _ => Polynomial { mons: mons, coefs: coefs }
    }
  }

  fn zero() -> Polynomial<M> {
    let one_mon: M = Monomial::one();
    Polynomial::new(~[one_mon], ~[0 as R])
  }

  fn domain_dim(&self) -> Dim {
    self.mons[0].domain_dim()
  }

  fn canonical_form(&self) -> Polynomial<M> {
    let mut coefs_by_mon: TreeMap<M,R> = TreeMap::new();
    for i in range(0, self.mons.len()) {
      let (mon, mon_coef) = (self.mons[i].clone(), self.coefs[i]);
      if mon_coef != 0 as R {
        let did_update = match coefs_by_mon.find_mut(&mon) {
          Some(coef) => { *coef += mon_coef; true }, None => false
        };
        if !did_update {
          coefs_by_mon.insert(mon, mon_coef);
        }
      }
    }
    let mut mons:  ~[M] = vec::with_capacity(coefs_by_mon.len());
    let mut coefs: ~[R] = vec::with_capacity(coefs_by_mon.len());
    for (mon, coef) in coefs_by_mon.iter() {
      if *coef != 0 as R {
        mons.push(mon.clone());
        coefs.push(*coef);
      }
    }
    if mons.len() != 0 {
      Polynomial::new(mons, coefs)
    }
    else {
      Polynomial::zero()
    }
  }
  
  fn equiv(&self, other: &Polynomial<M>) -> bool {
    let (can_self, can_other) = (self.canonical_form(), other.canonical_form());
    can_self.coefs == can_other.coefs && can_self.mons == can_other.mons
  }
  
  fn scale(&mut self, r: R) -> () {
    for i in range(0, self.coefs.len()) {
      self.coefs[i] *= r;
    }
  }
  
  fn scaled(&self, r: R) -> Polynomial<M> {
    let scaled_coefs = vec::from_fn(self.coefs.len(), |i| r * self.coefs[i]);
    Polynomial { mons: self.mons.clone(), coefs: scaled_coefs }
  }

}

impl<M:Monomial> Add<Polynomial<M>, Polynomial<M>> for Polynomial<M> {
  fn add(&self, other: &Polynomial<M>) -> Polynomial<M> {
    Polynomial {mons: self.mons + other.mons, coefs: self.coefs + other.coefs }
  }
}
impl<M:Monomial> Sub<Polynomial<M>, Polynomial<M>> for Polynomial<M> {
  fn sub(&self, other: &Polynomial<M>) -> Polynomial<M> {
    Polynomial {mons: self.mons + other.mons, coefs: self.coefs + other.coefs.map(|&c| -c)}
  }
}

impl<M:Monomial> Mul<Polynomial<M>, Polynomial<M>> for Polynomial<M> {
  fn mul(&self, other: &Polynomial<M>) -> Polynomial<M> {
    let (n1, n2) = (self.mons.len(), other.mons.len());
    let n = n1 * n2;
    let mut mons = vec::with_capacity(n);
    let mut coefs = vec::with_capacity(n);
    for i1 in range(0, n1) {
      for i2 in range(0, n2) {
        mons.push(self.mons[i1] * other.mons[i2]);
        coefs.push(self.coefs[i1] * other.coefs[i2]);
      }
    }
    Polynomial { mons: mons, coefs: coefs }
  }
}


impl<M:Monomial> ToStr for Polynomial<M> {
  fn to_str(&self) -> ~str {
    // join(map(i->"$(p.coefs[i]) $(p.mons[i])", 1:length(p.mons)), " + ")
    fmt!("Polynomial %?", *self)
  }
}


#[test]
fn test_1d() {
  // construction
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let x2_mon = Mon1d { exps: [Deg(2)] };
  let one = Polynomial::new(~[one_mon], ~[1.]);
  let zero: Polynomial<Mon1d> = Polynomial::zero();
  let zero2 = Polynomial::new(~[one_mon], ~[0.]);
  let x = Polynomial::new(~[x_mon], ~[1.]);
  let x2 = Polynomial::new(~[x_mon], ~[2.]);
  let one_plus_x = Polynomial::new(~[one_mon,x_mon], ~[1.,1.]);
  let two_plus_x = Polynomial::new(~[one_mon,x_mon], ~[2.,1.]);
  assert!(zero.equiv(&zero2));
  assert!(one_plus_x.equiv(&one_plus_x));
  assert!(!one_plus_x.equiv(&two_plus_x));

  // scaling
  assert!(one_plus_x.scaled(2.).equiv(&Polynomial::new(~[one_mon,x_mon], ~[2.,2.])));
  assert!(!one_plus_x.scaled(2.).equiv(&Polynomial::new(~[one_mon,x_mon], ~[2.,1.])));
  let mut p = one_plus_x.clone();
  p.scale(3.);
  assert!(p.equiv(&Polynomial::new(~[one_mon,x_mon], ~[3.,3.])));

  // addition and multiplication
  p = Polynomial::new(~[one_mon,x_mon,x2_mon], ~[1.,2.,1.]);
  assert!(p.equiv(&(one + x.scaled(2.) + x*x)));
  let sq_one_plus_x = one_plus_x * one_plus_x;
  assert!(p.equiv(&sq_one_plus_x))
}

/*
#[test]
fn test_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one: Mon2d = Monomial::one();
  let one_plus_x = Polynomial::new(~[one,x], ~[1 as R,1 as R]);
  let another_one_plus_x = Polynomial::new(~[one,x], ~[1 as R,1 as R]);
  let two_plus_x = Polynomial::new(~[one,x], ~[2 as R,1 as R]);
  // Test equality.
  assert_eq!(&one_plus_x, &one_plus_x);
  assert_eq!(&one_plus_x, &another_one_plus_x);
  assert!(one_plus_x == another_one_plus_x);
  assert!(one_plus_x != two_plus_x);
  assert!(x == x);
  assert!(x != y);
  
}
*/

/*
fn main() {
  let x = @Monomial::with_exps(~[deg(1), deg(0)]);
  let y = @Monomial::with_exps(~[deg(0), deg(1)]);
  let xy = @Monomial::with_exps(~[deg(1), deg(1)]);
  let mons = ~[x, y, xy, x, xy];
  let p1 = Polynomial::new(mons, ~[1., -1., 2., -1., 2.]);
  let p2 = Polynomial::new(~[xy, y], ~[1., -1.]);
  //let cfp = p1.canonical_form();
  //println(fmt!("Canonical form of %s is %s", p1.to_str(), cfp.to_str()));
  println(fmt!("p1 + p2 = %s, canon. form = %s", (p1 + p2).to_str(),(p1 + p2).canonical_form().to_str()));
  println(fmt!("p1 - p2 = %s, canon. form = %s", (p1 - p2).to_str(),(p1 - p2).canonical_form().to_str()));
  println(fmt!("2*(p1 - p2) =  %s", ((p1 - p2)*2.).canonical_form().to_str()));
}
*/

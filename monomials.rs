use common::*;
mod common;

pub trait Monomial: Eq +
                    TotalEq +
                    Ord +
                    TotalOrd +
                    IterBytes +
                    Clone +
                    ToStr +
                    Mul<Self,Self> {
                    // Mul<R, Self> { TODO: error "duplicate supertrait in trait def"

  fn domain_dim(&self) -> Dim;

  // fn exps_iter<'a>(&'a self) -> std::vec::VecIterator<'a, Deg>;

  fn map_exp(&self, n: Dim, f: &fn(Deg) -> Deg) -> Self;

  fn one() -> Self;
}


pub struct Mon1d {
  exps: [Deg,..1]
}
pub struct Mon2d {
  exps: [Deg,..2]
}
pub struct Mon3d {
  exps: [Deg,..3]
}
pub struct Mon4d {
  exps: [Deg,..4]
}

impl Monomial for Mon1d {
  fn domain_dim(&self) -> Dim {
    Dim(1)
  }
  /*
  fn exps_iter<'a>(&'a self) -> vec::VecIterator<'a,Deg> {
    self.exps.iter()
  }
  */
  fn map_exp(&self, n: Dim, f: &fn(Deg) -> Deg) -> Mon1d {
    if *n > 1 { fail!("Dimension number out of range in map_exp."); }
    else {
      Mon1d { exps: [ f(self.exps[0]) ] }
    }
  }

  fn one() -> Mon1d {
    Mon1d { exps: [Deg(0)] }
  }
}
impl Monomial for Mon2d {
  fn domain_dim(&self) -> Dim {
    Dim(2)
  }
  fn map_exp(&self, n: Dim, f: &fn(Deg) -> Deg) -> Mon2d {
    match n {
      Dim(0) => Mon2d { exps: [f(self.exps[0]),  self.exps[1] ] },
      Dim(1) => Mon2d { exps: [  self.exps[0], f(self.exps[1])] },
      _ => fail!("Dimension number out of range in map_exp.") 
    }
  }
  fn one() -> Mon2d {
    Mon2d { exps: [Deg(0), Deg(0)] }
  }
}
impl Monomial for Mon3d {
  fn domain_dim(&self) -> Dim {
    Dim(3)
  }
  fn map_exp(&self, n: Dim, f: &fn(Deg) -> Deg) -> Mon3d {
    match n {
      Dim(0) => Mon3d { exps: [f(self.exps[0]),  self.exps[1],   self.exps[2] ] },
      Dim(1) => Mon3d { exps: [  self.exps[0], f(self.exps[1]),  self.exps[2] ] },
      Dim(2) => Mon3d { exps: [  self.exps[0],   self.exps[1], f(self.exps[2])] },
      _ => fail!("Dimension number out of range in map_exp.") 
    }
  }
  fn one() -> Mon3d {
    Mon3d { exps: [Deg(0), Deg(0), Deg(0)] }
  }
}
impl Monomial for Mon4d {
  fn domain_dim(&self) -> Dim {
    Dim(4)
  }
  fn map_exp(&self, n: Dim, f: &fn(Deg) -> Deg) -> Mon4d {
    match n {
      Dim(0) => Mon4d { exps: [f(self.exps[0]),  self.exps[1],   self.exps[2],   self.exps[3] ] },
      Dim(1) => Mon4d { exps: [  self.exps[0], f(self.exps[1]),  self.exps[2],   self.exps[3] ] },
      Dim(2) => Mon4d { exps: [  self.exps[0],   self.exps[1], f(self.exps[2]),  self.exps[3] ] },
      Dim(3) => Mon4d { exps: [  self.exps[0],   self.exps[1],   self.exps[2], f(self.exps[3])] },
      _ => fail!("Dimension number out of range in map_exp.") 
    }
  }
  fn one() -> Mon4d {
    Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] }
  }
}

impl ToStr for Mon1d {
  fn to_str(&self) -> ~str {
    fmt!("x^%?", *self.exps[0])
  }
}
impl ToStr for Mon2d {
  fn to_str(&self) -> ~str {
    fmt!("x^%?y^%?", *self.exps[0], *self.exps[1])
  }
}
impl ToStr for Mon3d {
  fn to_str(&self) -> ~str {
    fmt!("x^%?y^%?z^%?", *self.exps[0], *self.exps[1], *self.exps[2])
  }
}
impl ToStr for Mon4d {
  fn to_str(&self) -> ~str {
    fmt!("x1^%?x2^%?x3^%?x4^%?", *self.exps[0], *self.exps[1], *self.exps[2], *self.exps[3])
  }
}

macro_rules! eq_impl(($t: ty) => {
  impl Eq for $t {
    fn eq(&self, other: &$t) -> bool {
      self.exps == other.exps
    }
    fn ne(&self, other: &$t) -> bool {
      self.exps != other.exps
    }
  }
})
eq_impl!(Mon1d)
eq_impl!(Mon2d)
eq_impl!(Mon3d)
eq_impl!(Mon4d)


// Implement IterBytes for hashing.
impl IterBytes for Mon1d {
  fn iter_bytes(&self, lsb0: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0]])
  }
}
impl IterBytes for Mon2d {
  fn iter_bytes(&self, lsb0: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0], *self.exps[1]])
  }
}
impl IterBytes for Mon3d {
  fn iter_bytes(&self, lsb0: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0], *self.exps[1], *self.exps[3]])
  }
}
impl IterBytes for Mon4d {
  fn iter_bytes(&self, lsb0: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0], *self.exps[1], *self.exps[3], *self.exps[4]])
  }
}

macro_rules! ord_impl(($t: ty) => {
  impl Ord for $t {
    fn lt(&self, other: &$t) -> bool {
      self.exps < other.exps
    }
    fn le(&self, other: &$t) -> bool {
      self.exps <= other.exps
    }
    fn gt(&self, other: &$t) -> bool {
      self.exps > other.exps
    }
    fn ge(&self, other: &$t) -> bool {
      self.exps >= other.exps
    }
  }
})
ord_impl!(Mon1d)
ord_impl!(Mon2d)
ord_impl!(Mon3d)
ord_impl!(Mon4d)

macro_rules! totord_impl(($t: ty) => {
  impl TotalOrd for $t {
    fn cmp(&self, other: &$t) -> Ordering {
      if self.exps < other.exps { Less }
      else if self.exps > other.exps { Greater }
      else { Equal }
    }
  }
})
totord_impl!(Mon1d)
totord_impl!(Mon2d)
totord_impl!(Mon3d)
totord_impl!(Mon4d)

macro_rules! toteq_impl(($t: ty) => {
  impl TotalEq for $t {
    fn equals(&self, other: &$t) -> bool {
      self.exps == other.exps
    }
  }
})
toteq_impl!(Mon1d)
toteq_impl!(Mon2d)
toteq_impl!(Mon3d)
toteq_impl!(Mon4d)

/* TODO: report macro bug
macro_rules! clone_impl(($t: ty) => {
  impl Clone for $t {
    fn clone(&self) -> $t {
      $t { exps: self.exps }
    }
  }
})
clone_impl!(Mon1d)
clone_impl!(Mon2d)
clone_impl!(Mon3d)
clone_impl!(Mon4d)
*/

impl Clone for Mon1d {
  fn clone(&self) -> Mon1d {
    Mon1d { exps: self.exps }
  }
}
impl Clone for Mon2d {
  fn clone(&self) -> Mon2d {
    Mon2d { exps: self.exps }
  }
}
impl Clone for Mon3d {
  fn clone(&self) -> Mon3d {
    Mon3d { exps: self.exps }
  }
} impl Clone for Mon4d {
  fn clone(&self) -> Mon4d {
    Mon4d { exps: self.exps }
  }
}


impl Mul<Mon1d, Mon1d> for Mon1d {
  fn mul(&self, other: &Mon1d) -> Mon1d {
    Mon1d { exps: [Deg(*self.exps[0] + *other.exps[0])] }
  }
}
impl Mul<Mon2d, Mon2d> for Mon2d {
  fn mul(&self, other: &Mon2d) -> Mon2d {
    Mon2d { exps: [Deg(*self.exps[0] + *other.exps[0]),
                   Deg(*self.exps[1] + *other.exps[1])] }
  }
}
impl Mul<Mon3d, Mon3d> for Mon3d {
  fn mul(&self, other: &Mon3d) -> Mon3d {
    Mon3d { exps: [Deg(*self.exps[0] + *other.exps[0]),
                   Deg(*self.exps[1] + *other.exps[1]),
                   Deg(*self.exps[2] + *other.exps[2])] }
  }
}
impl Mul<Mon4d, Mon4d> for Mon4d {
  fn mul(&self, other: &Mon4d) -> Mon4d {
    Mon4d { exps: [Deg(*self.exps[0] + *other.exps[0]),
                   Deg(*self.exps[1] + *other.exps[1]),
                   Deg(*self.exps[2] + *other.exps[2]),
                   Deg(*self.exps[3] + *other.exps[3])] }
  }
}

//////////////////////////////////////////////////////////
// Tests

#[test]
fn test_equality_1d() {
  let a1 = Mon1d { exps: [Deg(1)] };
  let b1 = Mon1d { exps: [Deg(1)] };
  assert_eq!(a1, b1);

  let c1 = Mon1d { exps: [Deg(0)] };
  assert!(a1 != c1);
  assert!(!(a1 == c1));
}

#[test]
fn test_equality_2d() {
  let a2 = Mon2d { exps: [Deg(1), Deg(0)] };
  let b2 = Mon2d { exps: [Deg(1), Deg(0)] };
  assert_eq!(a2, b2);

  let c2 = Mon2d { exps: [Deg(1), Deg(1)] };
  assert!(a2 != c2);
  assert!(!(a2 == c2));
}

#[test]
fn test_ord_1d() {
  let x = Mon1d { exps: [Deg(1)] };
  let one: Mon1d = Monomial::one();
  assert!(one < x);
  assert!(!(x < one));
  assert!(x > one);

  let x2 = Mon1d { exps: [Deg(2)] };
  assert!(x2 > x);
  assert!(x < x2);
  assert!(!(x > x2));
  assert!(!(x2 < x));
}

#[test]
fn test_ord_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let one: Mon2d = Monomial::one();
  assert!(one < x);
  assert!(!(x < one));
  assert!(x > one);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  let x2 = Mon2d { exps: [Deg(2), Deg(0)] };
  assert!(x2 > x);
  assert!(x < x2);
  assert!(!(x > x2));
  assert!(!(x2 < x));

  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  assert!(y < x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);
}

#[test]
fn test_clone_1d() {
  let x = Mon1d { exps: [Deg(1)] };
  assert_eq!(x, x.clone());
  let one = Mon1d { exps: [Deg(0)] };
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}

#[test]
fn test_clone_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  assert_eq!(x, x.clone());
  let one: Mon2d = Monomial::one();
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}
 
#[test]
fn test_hash_1d() {
  use std::hashmap::HashSet;
  let x = Mon1d { exps: [Deg(1)] };
  let one: Mon1d = Monomial::one();
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  m.insert(one);
  assert!(m.contains(&one));
}

#[test]
fn test_hash_2d() {
  use std::hashmap::HashSet;
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one: Mon2d = Monomial::one();
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  assert!(!m.contains(&y));
  m.insert(one);
  assert!(m.contains(&one));
  m.insert(y);
  assert!(m.contains(&y));
}

#[test]
fn test_mul_2d() {
  use std::hashmap::HashSet;
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one: Mon2d = Monomial::one();
  assert_eq!(x * y, Mon2d { exps: [Deg(1), Deg(1)] });
  assert_eq!(x * one, x);
  assert_eq!(one * y, y);
}

#[test]
fn test_mul_3d() {
  use std::hashmap::HashSet;
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one: Mon3d = Monomial::one();
  assert_eq!(x * y, Mon3d { exps: [Deg(1), Deg(1), Deg(0)] });
  assert_eq!(x * one, x);
  assert_eq!(one * y, y);
  assert_eq!(x*z, Mon3d { exps: [Deg(1), Deg(0), Deg(1)] });
  assert_eq!(x*y, Mon3d { exps: [Deg(1), Deg(1), Deg(0)] });
  assert_eq!(y*z, Mon3d { exps: [Deg(0), Deg(1), Deg(1)] });
}

#[test]
fn test_mul_4d() {
  use std::hashmap::HashSet;
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one: Mon4d = Monomial::one();
  assert_eq!(x * y, Mon4d { exps: [Deg(1), Deg(1), Deg(0), Deg(0)] });
  assert_eq!(x * one, x);
  assert_eq!(one * y, y);
  assert_eq!(x*z, Mon4d { exps: [Deg(1), Deg(0), Deg(1), Deg(0)] });
  assert_eq!(x*y, Mon4d { exps: [Deg(1), Deg(1), Deg(0), Deg(0)] });
  assert_eq!(y*z, Mon4d { exps: [Deg(0), Deg(1), Deg(1), Deg(0)] });

  assert_eq!(t*x, Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(1)] });
  assert_eq!(t*y, Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(1)] });
  assert_eq!(t*z, Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(1)] });
}

#[test]
fn test_map_exp_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one: Mon2d = Monomial::one();
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), x*x*y);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
}
#[test]
fn test_map_exp_3d() {
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one: Mon3d = Monomial::one();
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), x*x*y);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
  assert_eq!(z.map_exp(Dim(0), |e| Deg(*e+2)), z*x*x);
  assert_eq!(z.map_exp(Dim(1), |e| Deg(*e+2)), z*y*y);
  assert_eq!(z.map_exp(Dim(2), |e| Deg(*e+2)), z*z*z);
}


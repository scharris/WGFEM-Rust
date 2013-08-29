use common::*;
mod common;

trait Monomial {
  fn domain_dim(&self) -> Dim;
  // fn exps_iter<'a>(&'a self) -> std::vec::VecIterator<'a, Deg>;
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

impl Mon1d {
  fn with_exps(exps: [Deg,..1]) -> Mon1d {
    Mon1d { exps: exps }
  }
}
impl Mon2d {
  fn with_exps(exps: [Deg,..2]) -> Mon2d {
    Mon2d { exps: exps }
  }
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
  fn one() -> Mon1d {
    Mon1d { exps: [Deg(0)] }
  }
}
impl Monomial for Mon2d {
  fn domain_dim(&self) -> Dim {
    Dim(2)
  }
  fn one() -> Mon2d {
    Mon2d { exps: [Deg(0),Deg(0)] }
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

impl Eq for Mon1d {
  fn eq(&self, other: &Mon1d) -> bool {
    self.exps == other.exps
  }
  fn ne(&self, other: &Mon1d) -> bool {
    self.exps != other.exps
  }
}
impl Eq for Mon2d {
  fn eq(&self, other: &Mon2d) -> bool {
    self.exps == other.exps
  }
  fn ne(&self, other: &Mon2d) -> bool {
    self.exps != other.exps
  }
}

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

impl Ord for Mon1d {
  fn lt(&self, other: &Mon1d) -> bool {
    self.exps < other.exps
  }
  fn le(&self, other: &Mon1d) -> bool {
    self.exps <= other.exps
  }
  fn gt(&self, other: &Mon1d) -> bool {
    self.exps > other.exps
  }
  fn ge(&self, other: &Mon1d) -> bool {
    self.exps >= other.exps
  }
}
impl Ord for Mon2d {
  fn lt(&self, other: &Mon2d) -> bool {
    self.exps < other.exps
  }
  fn le(&self, other: &Mon2d) -> bool {
    self.exps <= other.exps
  }
  fn gt(&self, other: &Mon2d) -> bool {
    self.exps > other.exps
  }
  fn ge(&self, other: &Mon2d) -> bool {
    self.exps >= other.exps
  }
}

impl TotalOrd for Mon1d {
  fn cmp(&self, other: &Mon1d) -> Ordering {
    self.exps[0].cmp(&other.exps[0])
  }
}
impl TotalOrd for Mon2d {
  fn cmp(&self, other: &Mon2d) -> Ordering {
    if self.exps < other.exps { Less }
    else if self.exps > other.exps { Greater }
    else { Equal }
  }
}

impl TotalEq for Mon1d {
  fn equals(&self, other: &Mon1d) -> bool {
    self.exps == other.exps
  }
}
impl TotalEq for Mon2d {
  fn equals(&self, other: &Mon2d) -> bool {
    self.exps == other.exps
  }
}

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

impl Mul<Mon1d,Mon1d> for Mon1d {
  fn mul(&self, other: &Mon1d) -> Mon1d {
    Mon1d { exps: [Deg(*self.exps[0] + *other.exps[0])] }
  }
}
impl Mul<Mon2d,Mon2d> for Mon2d {
  fn mul(&self, other: &Mon2d) -> Mon2d {
    Mon2d { exps: [Deg(*self.exps[0] + *other.exps[0]),
                   Deg(*self.exps[1] + *other.exps[1])] }
  }
}

#[test]
fn test_equality_1d() {
    let a1 = Mon1d::with_exps([Deg(1)]);
    let b1 = Mon1d::with_exps([Deg(1)]);
    assert_eq!(a1, b1);

    let c1 = Mon1d::with_exps([Deg(0)]);
    assert!(a1 != c1);
    assert!(!(a1 == c1));
}

#[test]
fn test_equality_2d() {
    let a2 = Mon2d::with_exps([Deg(1),Deg(0)]);
    let b2 = Mon2d::with_exps([Deg(1),Deg(0)]);
    assert_eq!(a2, b2);

    let c2 = Mon2d::with_exps([Deg(1),Deg(1)]);
    assert!(a2 != c2);
    assert!(!(a2 == c2));
}

#[test]
fn test_ord_1d() {
    let x = Mon1d::with_exps([Deg(1)]);
    let one = Mon1d::with_exps([Deg(0)]);
    assert!(one < x);
    assert!(!(x < one));
    assert!(x > one);

    let x2 = Mon1d::with_exps([Deg(2)]);
    assert!(x2 > x);
    assert!(x < x2);
    assert!(!(x > x2));
    assert!(!(x2 < x));
}

#[test]
fn test_ord_2d() {
    let x = Mon2d::with_exps([Deg(1),Deg(0)]);
    let one = Mon2d::with_exps([Deg(0),Deg(0)]);
    assert!(one < x);
    assert!(!(x < one));
    assert!(x > one);

    assert!(one.cmp(&x) == Less);
    assert!(x.cmp(&one) == Greater);
    assert!(x.cmp(&x) == Equal);

    let x2 = Mon2d::with_exps([Deg(2),Deg(0)]);
    assert!(x2 > x);
    assert!(x < x2);
    assert!(!(x > x2));
    assert!(!(x2 < x));

    let y = Mon2d::with_exps([Deg(0),Deg(1)]);
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
    let x = Mon1d::with_exps([Deg(1)]);
    assert_eq!(x, x.clone());
    let one = Mon1d::with_exps([Deg(0)]);
    assert_eq!(one, one.clone());
    assert!(!(x.clone() == one));
}

#[test]
fn test_clone_2d() {
    let x = Mon2d::with_exps([Deg(1), Deg(0)]);
    assert_eq!(x, x.clone());
    let one = Mon2d::with_exps([Deg(0), Deg(0)]);
    assert_eq!(one, one.clone());
    assert!(!(x.clone() == one));
}
 
#[test]
fn test_hash_1d() {
  use std::hashmap::HashSet;
  let x = Mon1d::with_exps([Deg(1)]);
  let one = Mon1d::with_exps([Deg(0)]);
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
  let x = Mon2d::with_exps([Deg(1), Deg(0)]);
  let y = Mon2d::with_exps([Deg(0), Deg(1)]);
  let one = Mon2d::with_exps([Deg(0), Deg(0)]);
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



fn main() {
  let x1d = Mon1d::with_exps([Deg(1)]);
  let x2d = Mon2d::with_exps([Deg(1),Deg(0)]);
  let y2d = Mon2d::with_exps([Deg(0),Deg(1)]);
  printfln!("%?, %?", x1d.to_str(), x2d.to_str());
}

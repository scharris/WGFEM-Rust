use common::*;
use std::option::*;

pub trait Monomial: Eq +
                    TotalEq +
                    Ord +
                    TotalOrd +
                    IterBytes +
                    Clone +
                    ToStr +
                    Mul<Self,Self> {

  fn domain_dim(witness: Option<Self>) -> uint;
  
  fn value_at(&self, x: &[R]) -> R;
  
  fn value_at_for_origin(&self, x: &[R], origin: &[R]) -> R;
  
  fn exp(&self, coord: Dim) -> Deg;

  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Self;

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

static one_1d: Mon1d = Mon1d { exps: [Deg(0),..1] };
static one_2d: Mon2d = Mon2d { exps: [Deg(0),..2] };
static one_3d: Mon3d = Mon3d { exps: [Deg(0),..3] };
static one_4d: Mon4d = Mon4d { exps: [Deg(0),..4] };



impl Monomial for Mon1d {

  #[inline(always)]
  fn domain_dim(witness: Option<Mon1d>) -> uint {
    1u
  }

  #[inline(always)]
  fn value_at(&self, x: &[R]) -> R {
    pow(x[0], *self.exps[0] as uint)
  }
  
  #[inline(always)]
  fn value_at_for_origin(&self, x: &[R], origin: &[R]) -> R {
    pow(x[0] - origin[0], *self.exps[0] as uint)
  }

  #[inline(always)]
  fn exp(&self, coord: Dim) -> Deg {
    self.exps[*coord]
  }
  
  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Mon1d {
    if *coord > 1 { fail!("Dimension number out of range in map_exp."); }
    else {
      Mon1d { exps: [ f(self.exps[0]) ] }
    }
  }

  #[inline]
  fn one() -> Mon1d {
    one_1d
  }
}

impl Monomial for Mon2d {

  #[inline(always)]
  fn domain_dim(witness: Option<Mon2d>) -> uint {
    2u
  }
  
  #[inline(always)]
  fn value_at(&self, x: &[R]) -> R {
    pow(x[0], *self.exps[0] as uint) *
    pow(x[1], *self.exps[1] as uint)
  }
  
  #[inline(always)]
  fn value_at_for_origin(&self, x: &[R], origin: &[R]) -> R {
    pow(x[0] - origin[0], *self.exps[0] as uint) *
    pow(x[1] - origin[1], *self.exps[1] as uint)
  }

  #[inline(always)]
  fn exp(&self, coord: Dim) -> Deg {
    self.exps[*coord]
  }
  
  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Mon2d {
    match coord {
      Dim(0) => Mon2d { exps: [f(self.exps[0]),  self.exps[1] ] },
      Dim(1) => Mon2d { exps: [  self.exps[0], f(self.exps[1])] },
      _ => fail!("Dimension number out of range in map_exp.") 
    }
  }

  #[inline(always)]
  fn one() -> Mon2d {
    one_2d
  }
}

impl Monomial for Mon3d {

  #[inline(always)]
  fn domain_dim(witness: Option<Mon3d>) -> uint {
    3u
  }
  
  #[inline(always)]
  fn value_at(&self, x: &[R]) -> R {
    pow(x[0], *self.exps[0] as uint) *
    pow(x[1], *self.exps[1] as uint) *
    pow(x[2], *self.exps[2] as uint)
  }
  
  #[inline(always)]
  fn value_at_for_origin(&self, x: &[R], origin: &[R]) -> R {
    pow(x[0] - origin[0], *self.exps[0] as uint) *
    pow(x[1] - origin[1], *self.exps[1] as uint) *
    pow(x[2] - origin[2], *self.exps[2] as uint)
  }

  #[inline(always)]
  fn exp(&self, coord: Dim) -> Deg {
    self.exps[*coord]
  }

  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Mon3d {
    match coord {
      Dim(0) => Mon3d { exps: [f(self.exps[0]),  self.exps[1],   self.exps[2] ] },
      Dim(1) => Mon3d { exps: [  self.exps[0], f(self.exps[1]),  self.exps[2] ] },
      Dim(2) => Mon3d { exps: [  self.exps[0],   self.exps[1], f(self.exps[2])] },
      _ => fail!("Dimension number out of range in map_exp.") 
    }
  }

  #[inline(always)]
  fn one() -> Mon3d {
    one_3d
  }
}

impl Monomial for Mon4d {

  #[inline(always)]
  fn domain_dim(witness: Option<Mon4d>) -> uint {
    4u
  }

  #[inline(always)]
  fn value_at(&self, x: &[R]) -> R {
    pow(x[0], *self.exps[0] as uint) *
    pow(x[1], *self.exps[1] as uint) *
    pow(x[2], *self.exps[2] as uint) *
    pow(x[3], *self.exps[3] as uint)
  }
  
  #[inline(always)]
  fn value_at_for_origin(&self, x: &[R], origin: &[R]) -> R {
    pow(x[0] - origin[0], *self.exps[0] as uint) *
    pow(x[1] - origin[1], *self.exps[1] as uint) *
    pow(x[2] - origin[2], *self.exps[2] as uint) *
    pow(x[3] - origin[3], *self.exps[3] as uint)
  }

  #[inline(always)]
  fn exp(&self, coord: Dim) -> Deg {
    self.exps[*coord]
  }

  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Mon4d {
    match coord {
      Dim(0) => Mon4d { exps: [f(self.exps[0]),  self.exps[1],   self.exps[2],   self.exps[3] ] },
      Dim(1) => Mon4d { exps: [  self.exps[0], f(self.exps[1]),  self.exps[2],   self.exps[3] ] },
      Dim(2) => Mon4d { exps: [  self.exps[0],   self.exps[1], f(self.exps[2]),  self.exps[3] ] },
      Dim(3) => Mon4d { exps: [  self.exps[0],   self.exps[1],   self.exps[2], f(self.exps[3])] },
      _ => fail!("Dimension number out of range in map_exp.") 
    }
  }

  #[inline(always)]
  fn one() -> Mon4d {
    one_4d
  }
}


impl ToStr for Mon1d {
  fn to_str(&self) -> ~str {
    fmt!("x^%d",
         *self.exps[0] as int)
  }
}

impl ToStr for Mon2d {
  fn to_str(&self) -> ~str {
    fmt!("x^%dy^%d",
         *self.exps[0] as int, *self.exps[1] as int)
  }
}

impl ToStr for Mon3d {
  fn to_str(&self) -> ~str {
    fmt!("x^%dy^%dz^%d",
         *self.exps[0] as int, *self.exps[1] as int, *self.exps[2] as int)
  }
}

impl ToStr for Mon4d {
  fn to_str(&self) -> ~str {
    fmt!("x1^%dx2^%dx3^%dx4^%d",
         *self.exps[0] as int, *self.exps[1] as int, *self.exps[2] as int, *self.exps[3] as int)
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
    f(&[*self.exps[0], *self.exps[1], *self.exps[2]])
  }
}

impl IterBytes for Mon4d {
  fn iter_bytes(&self, lsb0: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0], *self.exps[1], *self.exps[2], *self.exps[3]])
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
}

impl Clone for Mon4d {
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


fn test_domain_dims() {
  assert_eq!(Monomial::domain_dim(None::<Mon1d>), 1u);
  assert_eq!(Monomial::domain_dim(None::<Mon2d>), 2u);
  assert_eq!(Monomial::domain_dim(None::<Mon3d>), 3u);
  assert_eq!(Monomial::domain_dim(None::<Mon4d>), 4u);
}

#[test]
fn test_value_at_1d() {
  let one = Mon1d { exps: [Deg(0)] };
  let x = Mon1d { exps: [Deg(1)] };
  let x2 = Mon1d { exps: [Deg(2)] };
  let x3 = Mon1d { exps: [Deg(3)] };
  assert_eq!(one.value_at(&[2.]), 1.);
  assert_eq!(x.value_at(&[2.]), 2.);
  assert_eq!(x2.value_at(&[2.]), 4.);
  assert_eq!(x3.value_at(&[2.]), 8.);
}
#[test]
fn test_value_at_2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  assert_eq!(one.value_at(&[2.,3.]), 1.);
  
  let x1y2 = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(x1y2.value_at(&[2.,3.]), 18.);
  
  let x3y1 = Mon2d { exps: [Deg(3), Deg(1)] };
  assert_eq!(x3y1.value_at(&[2.,3.]), 24.);
}
#[test]
fn test_value_at_3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at(&[2.,3.,4.]), 1.);
  
  let x1y2z3 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(x1y2z3.value_at(&[2.,3.,4.]), 18.*64.);
  
  let x3y1z2 = Mon3d { exps: [Deg(3), Deg(1), Deg(2)] };
  assert_eq!(x3y1z2.value_at(&[2.,3.,4.]), 24.*16.);
}
#[test]
fn test_value_at_4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at(&[2.,3.,4.,5.]), 1.);
  
  let x1y2z3t4 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(x1y2z3t4.value_at(&[2.,3.,4.,5.]), 18.*64.*5.*5.*5.*5.);
  
  let x3y1z2t4 = Mon4d { exps: [Deg(3), Deg(1), Deg(2), Deg(4)] };
  assert_eq!(x3y1z2t4.value_at(&[2.,3.,4.,5.]), 24.*16.*5.*5.*5.*5.);
}

#[test]
fn test_value_at_for_origin_1d() {
  let one = Mon1d { exps: [Deg(0)] };
  let x = Mon1d { exps: [Deg(1)] };
  let x2 = Mon1d { exps: [Deg(2)] };
  let x3 = Mon1d { exps: [Deg(3)] };

  assert_eq!(one.value_at_for_origin(&[3.],&[1.]), 1.);
  assert_eq!(x.value_at_for_origin(&[3.],&[1.]), 2.);
  assert_eq!(x2.value_at_for_origin(&[3.],&[1.]), 4.);
  assert_eq!(x3.value_at_for_origin(&[3.],&[1.]), 8.);
}
#[test]
fn test_value_at_for_origin_2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  assert_eq!(one.value_at_for_origin(&[2.,3.],&[12.,13.]), 1.);
  
  let x1y2 = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(x1y2.value_at_for_origin(&[2.,3.],&[-1.,-2.]), 3. * 25.);
  
  let x3y1 = Mon2d { exps: [Deg(3), Deg(1)] };
  assert_eq!(x3y1.value_at_for_origin(&[2.,3.],&[-1.,1.]), 27. * 2.);
}
#[test]
fn test_value_at_for_origin_3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at_for_origin(&[2.,3.,4.],&[1.,0.,2.]), 1.);
  
  let x1y2z3 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(x1y2z3.value_at_for_origin(&[2.,3.,4.],&[3.,1.,2.]), -1. * 4. * 8.);
  
  let x3y1z2 = Mon3d { exps: [Deg(3), Deg(1), Deg(2)] };
  assert_eq!(x3y1z2.value_at_for_origin(&[2.,3.,4.],&[3.,1.,2.]), -1. * 2. * 4.);
}
#[test]
fn test_value_at_for_origin_4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at_for_origin(&[2.,3.,4.,5.],&[1.,20.,30.,40.]), 1.);
  
  let x1y2z3t4 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(x1y2z3t4.value_at_for_origin(&[2.,3.,4.,5.],&[0.,5.,1.,2.]), 2. * 4. * 27. * 3.*3.*3.*3.);
  
  let x3y1z2t4 = Mon4d { exps: [Deg(3), Deg(1), Deg(2), Deg(4)] };
  assert_eq!(x3y1z2t4.value_at_for_origin(&[2.,3.,4.,5.],&[0.,5.,1.,2.]), 8. * (-2.) * 9. * 3.*3.*3.*3.);
}


#[test]
fn test_exp_1d() {
  let m = Mon1d { exps: [Deg(1)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
}

#[test]
#[should_fail]
fn test_exp_fail_1d() {
  let one: Mon1d = Monomial::one();
  let bad = one.exp(Dim(1)); 
}

#[test]
fn test_exp_2d() {
  let m = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
  assert_eq!(m.exp(Dim(1)), Deg(2));
}
#[test]
#[should_fail]
fn test_exp_fail_2d() {
  let one: Mon2d = Monomial::one();
  let bad = one.exp(Dim(2)); 
}
#[test]
fn test_exp_3d() {
  let m = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
  assert_eq!(m.exp(Dim(1)), Deg(2));
  assert_eq!(m.exp(Dim(2)), Deg(3));
}
#[test]
#[should_fail]
fn test_exp_fail_3d() {
  let one: Mon3d = Monomial::one();
  let bad = one.exp(Dim(3)); 
}
#[test]
fn test_exp_4d() {
  let m = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
  assert_eq!(m.exp(Dim(1)), Deg(2));
  assert_eq!(m.exp(Dim(2)), Deg(3));
  assert_eq!(m.exp(Dim(3)), Deg(4));
}
#[test]
#[should_fail]
fn test_exp_fail_4d() {
  let one: Mon4d = Monomial::one();
  let bad = one.exp(Dim(4)); 
}

#[test]
fn test_equality_1d() {
  let m1 = Mon1d { exps: [Deg(1)] };
  let m2 = Mon1d { exps: [Deg(1)] };
  assert_eq!(m1, m2);

  let m3 = Mon1d { exps: [Deg(0)] };
  assert!(m1 != m3);
  assert!(!(m1 == m3));
}
#[test]
fn test_equality_2d() {
  let m1 = Mon2d { exps: [Deg(1), Deg(2)] };
  let m2 = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(m1, m2);

  let m3 = Mon2d { exps: [Deg(1), Deg(0)] };
  let m4 = Mon2d { exps: [Deg(0), Deg(2)] };
  assert!(m1 != m3);
  assert!(m1 != m4);
  assert!(!(m1 == m3));
}
#[test]
fn test_equality_3d() {
  let m1 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  let m2 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(m1, m2);

  let m3 = Mon3d { exps: [Deg(0), Deg(2), Deg(3)] };
  let m4 = Mon3d { exps: [Deg(1), Deg(0), Deg(3)] };
  let m5 = Mon3d { exps: [Deg(1), Deg(2), Deg(0)] };
  assert!(m1 != m3);
  assert!(m1 != m4);
  assert!(m1 != m5);
  assert!(!(m1 == m3));
}
#[test]
fn test_equality_4d() {
  let m1 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  let m2 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(m1, m2);

  let m3 = Mon4d { exps: [Deg(0), Deg(2), Deg(3), Deg(4)] };
  let m4 = Mon4d { exps: [Deg(1), Deg(0), Deg(3), Deg(4)] };
  let m5 = Mon4d { exps: [Deg(1), Deg(2), Deg(0), Deg(4)] };
  let m6 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(0)] };
  assert!(m1 != m3);
  assert!(m1 != m4);
  assert!(m1 != m5);
  assert!(m1 != m6);
  assert!(!(m1 == m3));
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
  let one: Mon2d = Monomial::one();
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  assert!(one < y);
  assert!(x > one);
  assert!(!(x < one));

  assert!(y < x);
  assert!(y < x*y);
  assert!(x < x*y);
  assert!(x*y > x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  assert!(x*x > x);
  assert!(x < x*x);
  assert!(!(x > x*x));
  assert!(!(x*x < x));

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);
}
#[test]
fn test_ord_3d() {
  let one: Mon3d = Monomial::one();
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  assert!(one < y);
  assert!(x > one);
  assert!(!(x < one));

  assert!(y < x);
  assert!(y < x*y);
  assert!(x < x*y);
  assert!(x*y > x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  assert!(x*x > x);
  assert!(x < x*x);
  assert!(!(x > x*x));
  assert!(!(x*x < x));

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);

  assert!(z*z < x*y);
  assert!(y*z < x);
  assert!(x*y*z > x);
  assert!(z <= z);
}
#[test]
fn test_ord_4d() {
  let one: Mon4d = Monomial::one();
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  assert!(one < y);
  assert!(x > one);
  assert!(!(x < one));

  assert!(y < x);
  assert!(y < x*y);
  assert!(x < x*y);
  assert!(x*y > x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  assert!(x*x > x);
  assert!(x < x*x);
  assert!(!(x > x*x));
  assert!(!(x*x < x));

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);

  assert!(z*z < x*y);
  assert!(y*z < x);
  assert!(x*y*z > x);
  assert!(z <= z);

  assert!(y*z*t*t < x*y*z);
  assert!(t < z);
  assert!(t < y);
  assert!(y*z*t > t);
  assert!(t <= t);
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
fn test_clone_3d() {
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  assert_eq!(x, x.clone());
  let one: Mon3d = Monomial::one();
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}
#[test]
fn test_clone_4d() {
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  assert_eq!(x, x.clone());
  let one: Mon4d = Monomial::one();
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
fn test_hash_3d() {
  use std::hashmap::HashSet;
  let one: Mon3d = Monomial::one();
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  assert!(!m.contains(&y));
  assert!(!m.contains(&z));
  m.insert(one);
  assert!(m.contains(&one));
  m.insert(y);
  assert!(m.contains(&y));
  m.insert(z);
  assert!(m.contains(&z));
}
#[test]
fn test_hash_4d() {
  use std::hashmap::HashSet;
  let one: Mon4d = Monomial::one();
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  assert!(!m.contains(&y));
  assert!(!m.contains(&z));
  assert!(!m.contains(&t));
  m.insert(one);
  assert!(m.contains(&one));
  m.insert(y);
  assert!(m.contains(&y));
  m.insert(z);
  assert!(m.contains(&z));
  m.insert(t);
  assert!(m.contains(&t));
}


#[test]
fn test_mul_1d() {
  use std::hashmap::HashSet;
  let one: Mon1d = Monomial::one();
  let x = Mon1d { exps: [Deg(1)] };
  assert_eq!(x * x, Mon1d { exps: [Deg(2)] });
  assert_eq!(x * one, x);
}
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
  assert_eq!(x * y * z * t, Mon4d { exps: [Deg(1), Deg(1), Deg(1), Deg(1)] });
  assert_eq!(x * one, x);
  assert_eq!(x * one, x);
  assert_eq!(z * one, z);
  assert_eq!(t * one, t);
  assert_eq!(x*z, Mon4d { exps: [Deg(1), Deg(0), Deg(1), Deg(0)] });
  assert_eq!(x*y, Mon4d { exps: [Deg(1), Deg(1), Deg(0), Deg(0)] });
  assert_eq!(y*z, Mon4d { exps: [Deg(0), Deg(1), Deg(1), Deg(0)] });
  assert_eq!(t*z, Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(1)] });
  assert_eq!(t*y, Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(1)] });
  assert_eq!(t*x, Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(1)] });
}


#[test]
fn test_map_exp_1d() {
  let x = Mon1d { exps: [Deg(1)] };
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
}
#[test]
fn test_map_exp_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
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
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(x.map_exp(Dim(2), |e| Deg(*e+1)), x*z);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), x*x*y);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
  assert_eq!(y.map_exp(Dim(2), |e| Deg(*e+2)), y*z*z);
  assert_eq!(z.map_exp(Dim(0), |e| Deg(*e+3)), x*x*x*z);
  assert_eq!(z.map_exp(Dim(1), |e| Deg(*e+3)), y*y*y*z);
  assert_eq!(z.map_exp(Dim(2), |e| Deg(*e+3)), z*z*z*z);
}
#[test]
fn test_map_exp_4d() {
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(x.map_exp(Dim(2), |e| Deg(*e+1)), x*z);
  assert_eq!(x.map_exp(Dim(3), |e| Deg(*e+1)), x*t);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), y*x*x);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
  assert_eq!(y.map_exp(Dim(2), |e| Deg(*e+2)), y*z*z);
  assert_eq!(y.map_exp(Dim(3), |e| Deg(*e+2)), y*t*t);
  assert_eq!(z.map_exp(Dim(0), |e| Deg(*e+3)), z*x*x*x);
  assert_eq!(z.map_exp(Dim(1), |e| Deg(*e+3)), z*y*y*y);
  assert_eq!(z.map_exp(Dim(2), |e| Deg(*e+3)), z*z*z*z);
  assert_eq!(z.map_exp(Dim(3), |e| Deg(*e+3)), z*t*t*t);
  assert_eq!(t.map_exp(Dim(0), |e| Deg(*e+1)), t*x);
  assert_eq!(t.map_exp(Dim(1), |e| Deg(*e+1)), t*y);
  assert_eq!(t.map_exp(Dim(2), |e| Deg(*e+1)), t*z);
  assert_eq!(t.map_exp(Dim(3), |e| Deg(*e+1)), t*t);
}


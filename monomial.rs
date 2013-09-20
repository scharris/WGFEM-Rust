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


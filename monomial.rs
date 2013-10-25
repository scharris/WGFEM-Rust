use common::*;
use std::option::*;
use std::vec;
use std::iter::{range_inclusive, AdditiveIterator};
use std::num::pow_with_uint;


pub trait Monomial: Eq +
                    TotalEq +
                    Ord +
                    TotalOrd +
                    IterBytes +
                    Clone +
                    ToStr +
                    Mul<Self,Self> {

  fn domain_space_dims(_: Option<Self>) -> uint;
  
  fn value_at(&self, x: &[R]) -> R;
  
  fn value_at_for_origin(&self, x: &[R], origin: &[R]) -> R;

  fn value_at_reduced_dim_by_fixing(&self, reduced_dim_pt: &[R], fixed_dim: Dim, fixed_dim_val: R) -> R;
  
  fn exp(&self, coord: Dim) -> Deg;

  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Self;

  fn one() -> Self;

  fn mons_with_deg_lim_asc(deg_lim: DegLim) -> ~[Self];
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

pub enum DegLim {
  MaxMonDeg(u8),
  MaxMonFactorDeg(u8)
}

impl Monomial for Mon1d {

  #[inline(always)]
  fn domain_space_dims(_: Option<Mon1d>) -> uint {
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
  fn value_at_reduced_dim_by_fixing(&self, _reduced_dim_pt: &[R], fixed_dim: Dim, fixed_dim_val: R) -> R {
    if fixed_dim != Dim(0) { fail!("Dimension number out of range in map_exp."); }
    pow(fixed_dim_val, *self.exps[0] as uint)
  }

  #[inline(always)]
  fn exp(&self, coord: Dim) -> Deg {
    self.exps[*coord]
  }
  
  fn map_exp(&self, coord: Dim, f: &fn(Deg) -> Deg) -> Mon1d {
    if *coord > 0 { fail!("Dimension number out of range in map_exp."); }
    else {
      Mon1d { exps: [ f(self.exps[0]) ] }
    }
  }

  #[inline]
  fn one() -> Mon1d {
    one_1d
  }

  fn mons_with_deg_lim_asc(deg_lim: DegLim) -> ~[Mon1d] {
    match deg_lim {
      MaxMonDeg(deg) => vec::from_fn(deg as uint + 1, |e| Mon1d { exps: [ Deg(e as u8) ] }),
      MaxMonFactorDeg(deg) =>  vec::from_fn(deg as uint + 1, |e| Mon1d { exps: [ Deg(e as u8) ] })
    }  
  }
}

impl Monomial for Mon2d {

  #[inline(always)]
  fn domain_space_dims(_: Option<Mon2d>) -> uint {
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
  fn value_at_reduced_dim_by_fixing(&self, reduced_dim_pt: &[R], fixed_dim: Dim, fixed_dim_val: R) -> R {
    match fixed_dim {
      Dim(0) => 
        pow(fixed_dim_val,     *self.exps[0] as uint) *
        pow(reduced_dim_pt[0], *self.exps[1] as uint),
      Dim(1) => 
        pow(reduced_dim_pt[0], *self.exps[0] as uint) *
        pow(fixed_dim_val,     *self.exps[1] as uint),
      _ => fail!("Dimension number out of range in map_exp.")
    }
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

  fn mons_with_deg_lim_asc(deg_lim: DegLim) -> ~[Mon2d] {
    match deg_lim {
      MaxMonDeg(deg) => {
        let mut mons: ~[Mon2d] = vec::with_capacity(num_mons_of_deg_le(Deg(deg), 2));
        for e0 in range_inclusive(0, deg) {
          for e1 in range_inclusive(0, deg - e0) {
            mons.push(Mon2d { exps: [Deg(e0), Deg(e1)] });
          }
        }
        mons
      }
      MaxMonFactorDeg(deg) => {
        let mut mons: ~[Mon2d] = vec::with_capacity(sq(deg as uint + 1));
        for e0 in range_inclusive(0, deg) {
          for e1 in range_inclusive(0, deg) {
            mons.push(Mon2d { exps: [Deg(e0), Deg(e1)] });
          }
        }
        mons
      }
    }  
  }

}

impl Monomial for Mon3d {

  #[inline(always)]
  fn domain_space_dims(_: Option<Mon3d>) -> uint {
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
  fn value_at_reduced_dim_by_fixing(&self, reduced_dim_pt: &[R], fixed_dim: Dim, fixed_dim_val: R) -> R {
    match fixed_dim {
      Dim(0) => 
        pow(fixed_dim_val,     *self.exps[0] as uint) *
        pow(reduced_dim_pt[0], *self.exps[1] as uint) *
        pow(reduced_dim_pt[1], *self.exps[2] as uint),
      Dim(1) => 
        pow(reduced_dim_pt[0], *self.exps[0] as uint) *
        pow(fixed_dim_val,     *self.exps[1] as uint) *
        pow(reduced_dim_pt[1], *self.exps[2] as uint),
      Dim(2) => 
        pow(reduced_dim_pt[0], *self.exps[0] as uint) *
        pow(reduced_dim_pt[1], *self.exps[1] as uint) *
        pow(fixed_dim_val,     *self.exps[2] as uint),
      _ => fail!("Dimension number out of range in map_exp.")
    }
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

  fn mons_with_deg_lim_asc(deg_lim: DegLim) -> ~[Mon3d] {
    match deg_lim {
      MaxMonDeg(deg) => {
        let mut mons: ~[Mon3d] = vec::with_capacity(num_mons_of_deg_le(Deg(deg), 3));
        for e0 in range_inclusive(0, deg) {
          for e1 in range_inclusive(0, deg - e0) {
            for e2 in range_inclusive(0, deg - e0 - e1) {
              mons.push(Mon3d { exps: [Deg(e0), Deg(e1), Deg(e2)] });
            }
          }
        }
        mons
      }
      MaxMonFactorDeg(deg) => {
        let mut mons: ~[Mon3d] = vec::with_capacity(pow_with_uint(deg as uint + 1, 3));
        for e0 in range_inclusive(0, deg) {
          for e1 in range_inclusive(0, deg) {
            for e2 in range_inclusive(0, deg) {
              mons.push(Mon3d { exps: [Deg(e0), Deg(e1), Deg(e2)] });
            }
          }
        }
        mons
      }
    }  
  }
}

impl Monomial for Mon4d {

  #[inline(always)]
  fn domain_space_dims(_: Option<Mon4d>) -> uint {
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
  fn value_at_reduced_dim_by_fixing(&self, reduced_dim_pt: &[R], fixed_dim: Dim, fixed_dim_val: R) -> R {
    match fixed_dim {
      Dim(0) => 
        pow(fixed_dim_val,     *self.exps[0] as uint) *
        pow(reduced_dim_pt[0], *self.exps[1] as uint) *
        pow(reduced_dim_pt[1], *self.exps[2] as uint) *
        pow(reduced_dim_pt[2], *self.exps[3] as uint),
      Dim(1) => 
        pow(reduced_dim_pt[0], *self.exps[0] as uint) *
        pow(fixed_dim_val,     *self.exps[1] as uint) *
        pow(reduced_dim_pt[1], *self.exps[2] as uint) *
        pow(reduced_dim_pt[2], *self.exps[3] as uint),
      Dim(2) => 
        pow(reduced_dim_pt[0], *self.exps[0] as uint) *
        pow(reduced_dim_pt[1], *self.exps[1] as uint) *
        pow(fixed_dim_val,     *self.exps[2] as uint) *
        pow(reduced_dim_pt[2], *self.exps[3] as uint),
      Dim(3) => 
        pow(reduced_dim_pt[0], *self.exps[0] as uint) *
        pow(reduced_dim_pt[1], *self.exps[1] as uint) *
        pow(reduced_dim_pt[2], *self.exps[2] as uint) *
        pow(fixed_dim_val,     *self.exps[3] as uint),
      _ => fail!("Dimension number out of range in map_exp.")
    }
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

  fn mons_with_deg_lim_asc(deg_lim: DegLim) -> ~[Mon4d] {
    match deg_lim {
      MaxMonDeg(deg) => {
        let mut mons: ~[Mon4d] = vec::with_capacity(num_mons_of_deg_le(Deg(deg), 4));
        for e0 in range_inclusive(0, deg) {
          for e1 in range_inclusive(0, deg - e0) {
            for e2 in range_inclusive(0, deg - e0 - e1) {
              for e3 in range_inclusive(0, deg - e0 - e1 - e2) {
                mons.push(Mon4d { exps: [Deg(e0), Deg(e1), Deg(e2), Deg(e3)] });
              }
            }
          }
        }
        mons
      }
      MaxMonFactorDeg(deg) => {
        let mut mons: ~[Mon4d] = vec::with_capacity(pow_with_uint(deg as uint + 1, 4));
        for e0 in range_inclusive(0, deg) {
          for e1 in range_inclusive(0, deg) {
            for e2 in range_inclusive(0, deg) {
              for e3 in range_inclusive(0, deg) {
                mons.push(Mon4d { exps: [Deg(e0), Deg(e1), Deg(e2), Deg(e3)] });
              }
            }
          }
        }
        mons
      }
    }  
  }
}


impl ToStr for Mon1d {
  fn to_str(&self) -> ~str {
    format!("x^{:d}",
            *self.exps[0] as int)
  }
}

impl ToStr for Mon2d {
  fn to_str(&self) -> ~str {
    format!("x^{:d}y^{:d}",
            *self.exps[0] as int, *self.exps[1] as int)
  }
}

impl ToStr for Mon3d {
  fn to_str(&self) -> ~str {
    format!("x^{:d}y^{:d}z^{:d}",
            *self.exps[0] as int, *self.exps[1] as int, *self.exps[2] as int)
  }
}

impl ToStr for Mon4d {
  fn to_str(&self) -> ~str {
    format!("x1^{:d}x2^{:d}x3^{:d}x4^{:d}",
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
  fn iter_bytes(&self, _: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0]])
  }
}

impl IterBytes for Mon2d {
  fn iter_bytes(&self, _: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0], *self.exps[1]])
  }
}

impl IterBytes for Mon3d {
  fn iter_bytes(&self, _: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
    f(&[*self.exps[0], *self.exps[1], *self.exps[2]])
  }
}

impl IterBytes for Mon4d {
  fn iter_bytes(&self, _: bool, f: &fn(buf: &[u8]) -> bool) -> bool {
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


static one_1d: Mon1d = Mon1d { exps: [Deg(0),..1] };
static one_2d: Mon2d = Mon2d { exps: [Deg(0),..2] };
static one_3d: Mon3d = Mon3d { exps: [Deg(0),..3] };
static one_4d: Mon4d = Mon4d { exps: [Deg(0),..4] };


// auxiliary functions


fn num_mons_of_deg_eq(deg: Deg, dom_space_dims: uint) -> uint {
  multiset_choose(dom_space_dims, *deg as uint)
}

pub fn num_mons_of_deg_le(max_deg: Deg, dom_space_dims: uint) -> uint {
  range_inclusive(0, *max_deg)
    .map(|deg| num_mons_of_deg_eq(Deg(deg), dom_space_dims))
    .sum()
}

pub fn num_mons_with_deg_lim(deg_lim: DegLim, dom_space_dims: uint) -> uint {
  match deg_lim {
    MaxMonDeg(max_deg) => range_inclusive(0, max_deg).map(|deg| num_mons_of_deg_eq(Deg(deg), dom_space_dims)).sum(),
    MaxMonFactorDeg(max_deg) => pow_with_uint(max_deg as uint + 1, dom_space_dims)
  }
}

#[inline]
pub fn domain_space_dims<Mon:Monomial>() -> uint {
  Monomial::domain_space_dims(None::<Mon>)
}

#[inline]
fn multiset_choose(n: uint, k: uint) -> uint {
  binomial(n + k - 1, k)
}

fn binomial(n: uint, k: uint) -> uint {
  if k == 0u || n == k { 1u }
  else if k > n { 0u }
  else {
    let k = if k > n - k { n - k } else { k };
    if k == 1 { n }
    else {
      let mut b = 1;
      for i in range_inclusive(1,k) {
        b *= (n - (k - i));
        b /= i;
      }
      b
    }
  }
}


#[test]
fn test_binomial() {
  assert_eq!(binomial(5, 0), 1);
  assert_eq!(binomial(5, 5), 1); 
  assert_eq!(binomial(5, 1), 5);
  assert_eq!(binomial(5, 2), 10);
}


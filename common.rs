use std::num::{Num,Zero};
use std::vec;
use std::f64;

// types and type aliases

pub type R = f64;
pub static R_NaN: R = f64::NAN;

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Deg(u8);

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Dim(uint);

// constants

pub static DEFAULT_INTEGRATION_REL_ERR: R = 1e-12;
pub static DEFAULT_INTEGRATION_ABS_ERR: R = 1e-12;


#[inline]
pub fn pow(radix: R, exp: uint) -> R {
  if exp == 0u { return 1 as R; }
  if radix == 0. { return 0 as R; }
  let mut prod = 1 as R;
  let mut _exp = exp;
  let mut multiplier = radix;
  while (_exp > 0u) {
    if _exp % 2u == 1u {
      prod *= multiplier;
    }
    _exp /= 2u;
    multiplier *= multiplier;
  }
 prod 
}

#[inline(always)]
pub fn sq<T:Mul<T,T>>(x: T) -> T {
  x.mul(&x)
}

pub fn cumulative_sums_prev_elems<T:Num+Clone>(xs: &[T]) -> ~[T] {
  let mut cum_sums = vec::with_capacity(xs.len());
  let mut sum: T =  Zero::zero();
  for x in xs.iter() {
    cum_sums.push(sum.clone());
    sum = sum.add(x);
  }
  cum_sums
}

#[test]
fn test_pow() {
  assert_eq!(pow(2.,2), 4.);
  assert_eq!(pow(2.,3), 8.);
  assert_eq!(pow(0.,3), 0.);
  assert_eq!(pow(1.,3), 1.);
  assert_eq!(pow(0.,0), 1.); // In context of integral exponents and continuous radix, this def is most appropriate.
}

#[test]
fn test_sq() {
  assert_eq!(sq(0.), 0.);
  assert_eq!(sq(1.), 1.);
  assert_eq!(sq(2.), 4.);
  assert_eq!(sq(2.2), 2.2*2.2);
}

#[test]
fn test_cum_sums_prev() {
  let xs =  ~[1u, 2u, 3u];
  assert_eq!(cumulative_sums_prev_elems(xs), ~[0u, 1u, 3u]);
}


// types and type aliases

pub type R = f64;

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Deg(u8);

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Dim(uint);

// constants

pub static DEFAULT_INTEGRATION_REL_ERR: R = 1e-12;
pub static DEFAULT_INTEGRATION_ABS_ERR: R = 1e-12;

// utility functions

#[inline(always)]
pub fn pow(base: R, exp: uint) -> R {
  let mut prod = 1 as R;
  let mut i = exp;
  while i != 0u {
    prod *= base;
    i -= 1;
  } 
  prod
}

#[test]
fn test_pow() {
  assert_eq!(pow(2f64, 0), 1f64);
  assert_eq!(pow(2f64, 1), 2f64);
  assert_eq!(pow(2f64, 2), 4f64);
}


// types and type aliases

pub type R = f64;

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

pub fn sq(x: uint) -> uint {
  x * x
}

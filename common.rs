
pub type R = f64;

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Deg(u8);

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Dim(uint);

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

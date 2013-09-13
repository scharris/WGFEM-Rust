
pub type R = f64;

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Deg(u8);

#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Dim(uint);

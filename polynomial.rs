use extra::treemap::TreeMap;
use std::vec;
use std::iter::Zip;
use common::*;
use monomial::*;

#[deriving(Eq, Clone)]
pub struct PolyWithOwnedMons<M> {
  coefs: ~[R],
  mons:  ~[M]
}

#[deriving(Eq, Clone)]
pub struct PolyWithBorrowedMons<'self,M> {
  coefs: ~[R],
  mons: &'self [M]
}

// Convenience function to create polynomials in testing code. Performance critical code
// should use the ::new function implementations instead which will be more efficient.
pub fn poly<M:Monomial>(terms: ~[(R,M)]) -> PolyWithOwnedMons<M> {
  let (coefs, mons) = vec::unzip(terms.move_iter());
  PolyWithOwnedMons { coefs: coefs, mons: mons } 
}

impl<M:Monomial> PolyWithOwnedMons<M> {
  
  pub fn new(coefs: ~[R], mons: ~[M]) -> PolyWithOwnedMons<M> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyWithOwnedMons { coefs: coefs, mons: mons }
    }
  }
  
  pub fn zero() -> PolyWithOwnedMons<M> {
    let one_mon: M = Monomial::one();
    PolyWithOwnedMons { coefs: ~[0 as R], mons: ~[one_mon] }
  }
  
  pub fn zero_with_capacity(n: uint) -> PolyWithOwnedMons<M> {
    let one_mon: M = Monomial::one();
    let mut coefs = vec::with_capacity(n);
    let mut mons = vec::with_capacity(n);
    coefs.push(0 as R);
    mons.push(one_mon);
    PolyWithOwnedMons { coefs: coefs, mons: mons }
  }
}


impl<'self,M:Monomial> PolyWithBorrowedMons<'self,M> {
  
  pub fn new(coefs: ~[R], mons: &'self [M]) -> PolyWithBorrowedMons<'self,M> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyWithBorrowedMons { coefs: coefs, mons: mons }
    }
  }
}

pub trait Polynomial<M>: Add<Self,PolyWithOwnedMons<M>> +
                         Mul<Self,PolyWithOwnedMons<M>> {

  fn domain_dim(witness: Option<Self>) -> uint;

  fn num_terms(&self) -> uint;
  
  fn term(&self, n: uint) -> (R,M);

  fn foldl_terms<A>(&self, z: A, f: &fn(a: A, term: (R,M)) -> A) -> A;
  
  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  ();

  fn scaled(&self, r: R) -> Self;

  fn scale(&mut self, r: R) -> ();

  fn canonical_form(&self) -> PolyWithOwnedMons<M>;
  
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool;

}


fn canonical_form_impl<M: Monomial>(coefs: &[R], mons: &[M]) -> (~[R], ~[M]) {
  let mut coefs_by_mon: TreeMap<M,R> = TreeMap::new();
  for i in range(0, mons.len()) {
    let (mon, mon_coef) = (mons[i].clone(), coefs[i]);
    if mon_coef != 0 as R {
      let did_update = match coefs_by_mon.find_mut(&mon) {
        Some(coef) => { *coef += mon_coef; true }, None => false
      };
      if !did_update {
        coefs_by_mon.insert(mon, mon_coef);
      }
    }
  }
  let mut mons:  ~[M] = vec::with_capacity(coefs_by_mon.len());
  let mut coefs: ~[R] = vec::with_capacity(coefs_by_mon.len());
  for (mon, coef) in coefs_by_mon.iter() {
    if *coef != 0 as R {
      mons.push(mon.clone());
      coefs.push(*coef);
    }
  }
  (coefs, mons)
}


impl<'self,M:Monomial> Polynomial<M>
                   for PolyWithBorrowedMons<'self,M> {

  #[inline]
  fn domain_dim(witness: Option<PolyWithBorrowedMons<M>>) -> uint {
    Monomial::domain_dim(None::<M>)
  }
  
  #[inline]
  fn num_terms(&self) -> uint {
    self.mons.len()
  }

  #[inline]
  fn term(&self, n: uint) -> (R,M) {
    let coef = self.coefs[n];
    let mon = unsafe { self.mons.unsafe_get(n) }; // bounds-checked already in coef access
    (coef,mon)
  }
  
  fn foldl_terms<A>(&self, z: A, f: &fn(a: A, term: (R,M)) -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.mons.len()) {
      let term = unsafe { (self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }

  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn scaled(&self, r: R) -> PolyWithBorrowedMons<'self, M> {
    let scaled_coefs = vec::from_fn(self.coefs.len(), |i| r * self.coefs[i]);
    PolyWithBorrowedMons { coefs: scaled_coefs, mons: self.mons }
  }
  
  fn scale(&mut self, r: R) -> () {
    for i in range(0, self.coefs.len()) {
      self.coefs[i] *= r;
    }
  }
  
  fn canonical_form(&self) -> PolyWithOwnedMons<M> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let zero: PolyWithOwnedMons<M> = PolyWithOwnedMons::zero(); zero}
      _ => PolyWithOwnedMons::new(can_coefs, can_mons)
    }
  }
  
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }
  
}

impl<M:Monomial> Polynomial<M>
             for PolyWithOwnedMons<M> {

  fn scaled(&self, r: R) -> PolyWithOwnedMons<M> {
    let scaled_coefs = vec::from_fn(self.coefs.len(), |i| r * self.coefs[i]);
    PolyWithOwnedMons { coefs: scaled_coefs, mons: self.mons.clone() }
  }
  
  fn scale(&mut self, r: R) -> () {
    for i in range(0, self.coefs.len()) {
      self.coefs[i] *= r;
    }
  }

  #[inline]
  fn domain_dim(witness: Option<PolyWithOwnedMons<M>>) -> uint {
    Monomial::domain_dim(None::<M>)
  }
  
  #[inline]
  fn num_terms(&self) -> uint {
    self.mons.len()
  }
  
  #[inline]
  fn term(&self, n: uint) -> (R,M) {
    let coef = self.coefs[n];
    let mon = unsafe { self.mons.unsafe_get(n) }; // bounds-checked already in coef access
    (coef,mon)
  }
  
  fn foldl_terms<A>(&self, z: A, f: &fn(a: A, term: (R,M)) -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.mons.len()) {
      // let term = (self.coefs[n], self.mons[n].clone());
      let term = unsafe { (self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }

  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyWithOwnedMons<M> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let m: PolyWithOwnedMons<M> = PolyWithOwnedMons::zero(); m}
      _ => PolyWithOwnedMons::new(can_coefs, can_mons)
    }
  }
 
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }

}

fn mul_polys_impl<M:Monomial>(coefs1: &[R], mons1: &[M], coefs2: &[R], mons2: &[M]) -> PolyWithOwnedMons<M> {
  let (n1, n2) = (mons1.len(), mons2.len());
  let n = n1 * n2;
  let mut mons = vec::with_capacity(n);
  let mut coefs = vec::with_capacity(n);
  for i1 in range(0, n1) {
    for i2 in range(0, n2) {
      mons.push(mons1[i1] * mons2[i2]);
      coefs.push(coefs1[i1] * coefs2[i2]);
    }
  }
  PolyWithOwnedMons { coefs: coefs, mons: mons }
}

impl<'self,M:Monomial> Mul<PolyWithBorrowedMons<'self,M>, PolyWithOwnedMons<M>>
                   for PolyWithBorrowedMons<'self,M> {

  #[inline]
  fn mul(&self, other: &PolyWithBorrowedMons<'self,M>) -> PolyWithOwnedMons<M> {
    mul_polys_impl(self.coefs, self.mons, other.coefs, other.mons)
  }
}

impl<M:Monomial> Mul<PolyWithOwnedMons<M>, PolyWithOwnedMons<M>>
             for PolyWithOwnedMons<M> {
    
  #[inline]
  fn mul(&self, other: &PolyWithOwnedMons<M>) -> PolyWithOwnedMons<M> {
    mul_polys_impl(self.coefs, self.mons, other.coefs, other.mons)
  }
}

impl<'self,M:Monomial> Add<PolyWithBorrowedMons<'self,M>, PolyWithOwnedMons<M>>
                   for PolyWithBorrowedMons<'self,M> {

  #[inline]
  fn add(&self, other: &PolyWithBorrowedMons<'self,M>) -> PolyWithOwnedMons<M> {
    PolyWithOwnedMons::new(self.coefs + other.coefs, self.mons + other.mons) 
  }
}

impl<M:Monomial> Add<PolyWithOwnedMons<M>, PolyWithOwnedMons<M>>
             for PolyWithOwnedMons<M> {
    
  #[inline]
  fn add(&self, other: &PolyWithOwnedMons<M>) -> PolyWithOwnedMons<M> {
    PolyWithOwnedMons::new(self.coefs + other.coefs, self.mons + other.mons)
  }
}

impl<'self,M:Monomial> Sub<PolyWithBorrowedMons<'self,M>, PolyWithOwnedMons<M>>
                   for PolyWithBorrowedMons<'self,M> {

  #[inline]
  fn sub(&self, other: &PolyWithBorrowedMons<'self,M>) -> PolyWithOwnedMons<M> {
    PolyWithOwnedMons::new(self.coefs + other.coefs.map(|&c| -c), self.mons + other.mons) 
  }
}

impl<M:Monomial> Sub<PolyWithOwnedMons<M>, PolyWithOwnedMons<M>>
             for PolyWithOwnedMons<M> {
    
  #[inline]
  fn sub(&self, other: &PolyWithOwnedMons<M>) -> PolyWithOwnedMons<M> {
    PolyWithOwnedMons::new(self.coefs + other.coefs.map(|&c| -c), self.mons + other.mons) 
  }
}

fn to_str_impl<M: Monomial>(coefs: &[R], mons: &[M]) -> ~str {
  let term_strs = vec::from_fn(mons.len(), |i| {
    coefs[i].to_str() + " " + mons[i].to_str()
  });
  term_strs.connect(" + ")
}


impl<'self,M:Monomial> ToStr
                   for PolyWithBorrowedMons<'self,M> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

impl<M:Monomial> ToStr
             for PolyWithOwnedMons<M> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}


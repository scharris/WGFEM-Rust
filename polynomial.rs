use extra::treemap::TreeMap;
use std::vec;
use std::num::abs;

use common::*;
use monomial::*;


#[deriving(Eq, Clone)]
pub struct PolyOwning<M> {
  coefs: ~[R],
  mons:  ~[M]
}

#[deriving(Eq, Clone)]
pub struct PolyBorrowing<'self,M> {
  coefs: &'self [R],
  mons: &'self [M]
}

#[deriving(Eq, Clone)]
pub struct PolyBorrowingMons<'self,M> {
  coefs: ~[R],
  mons: &'self [M]
}


impl<M:Monomial> PolyOwning<M> {
  
  pub fn new(coefs: ~[R], mons: ~[M]) -> PolyOwning<M> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyOwning { coefs: coefs, mons: mons }
    }
  }
  
  pub fn zero() -> PolyOwning<M> {
    let one_mon: M = Monomial::one();
    PolyOwning { coefs: ~[0 as R], mons: ~[one_mon] }
  }
  
  pub fn zero_with_capacity(n: uint) -> PolyOwning<M> {
    let one_mon: M = Monomial::one();
    let mut coefs = vec::with_capacity(n);
    let mut mons = vec::with_capacity(n);
    coefs.push(0 as R);
    mons.push(one_mon);
    PolyOwning { coefs: coefs, mons: mons }
  }
  
  pub fn from_polys_lcomb<P:Polynomial<M>>(terms: &[(R,&P)]) -> PolyOwning<M> {
    let mut coefs_by_mon: TreeMap<M,R> = TreeMap::new();
    for &(c_p, p) in terms.iter() {
      p.each_term(|(c_m, m)| {
        if c_m != 0 as R {
          let did_update = match coefs_by_mon.find_mut(&m) {
            Some(c) => { *c += c_p * c_m; true }, None => false
          };
          if !did_update {
            coefs_by_mon.insert(m, c_p * c_m);
          }
        }
      });
    }
    let mut mons:  ~[M] = vec::with_capacity(coefs_by_mon.len());
    let mut coefs: ~[R] = vec::with_capacity(coefs_by_mon.len());
    for (mon, coef) in coefs_by_mon.iter() {
      if *coef != 0 as R {
        mons.push(mon.clone());
        coefs.push(*coef);
      }
    }
    PolyOwning { coefs: coefs, mons: mons }
  }

}

impl<'self,M:Monomial> PolyBorrowing<'self,M> {
  
  pub fn new(coefs: &'self [R], mons: &'self [M]) -> PolyBorrowing<'self,M> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyBorrowing { coefs: coefs, mons: mons }
    }
  }
}

impl<'self,M:Monomial> PolyBorrowingMons<'self,M> {
  
  pub fn new(coefs: ~[R], mons: &'self [M]) -> PolyBorrowingMons<'self,M> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyBorrowingMons { coefs: coefs, mons: mons }
    }
  }
}

pub trait Polynomial<M>: ToStr {

  fn domain_space_dims(_: Option<Self>) -> uint;

  fn num_terms(&self) -> uint;
  
  fn term(&self, n: uint) -> (R,M);

  fn foldl_terms<A>(&self, z: A, f: &fn(a: A, term: (R,M)) -> A) -> A;
  
  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  ();

  fn canonical_form(&self) -> PolyOwning<M>;
  
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool;

}


fn canonical_form_impl<M: Monomial>(coefs: &[R], mons: &[M]) -> (~[R],~[M]) {
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


impl<M:Monomial> Polynomial<M>
             for PolyOwning<M> {

  #[inline]
  fn domain_space_dims(_: Option<PolyOwning<M>>) -> uint {
    Monomial::domain_space_dims(None::<M>)
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

  #[inline]
  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyOwning<M> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let m: PolyOwning<M> = PolyOwning::zero(); m}
      _ => PolyOwning::new(can_coefs, can_mons)
    }
  }
 
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }

}

impl<'self,M:Monomial> Polynomial<M>
                   for PolyBorrowing<'self,M> {

  #[inline]
  fn domain_space_dims(_: Option<PolyBorrowing<M>>) -> uint {
    Monomial::domain_space_dims(None::<M>)
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

  #[inline]
  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyOwning<M> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let zero: PolyOwning<M> = PolyOwning::zero(); zero}
      _ => PolyOwning::new(can_coefs, can_mons)
    }
  }
  
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }
 
}

impl<'self,M:Monomial> Polynomial<M>
                   for PolyBorrowingMons<'self,M> {

  #[inline]
  fn domain_space_dims(_: Option<PolyBorrowingMons<M>>) -> uint {
    Monomial::domain_space_dims(None::<M>)
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

  #[inline]
  fn each_term(&self, f: &fn(term: (R,M)) -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyOwning<M> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let zero: PolyOwning<M> = PolyOwning::zero(); zero}
      _ => PolyOwning::new(can_coefs, can_mons)
    }
  }
  
  fn equiv<P: Polynomial<M>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }
 
}

// Implement ToStr trait.

fn to_str_impl<M: Monomial>(coefs: &[R], mons: &[M]) -> ~str {
  let term_strs = vec::from_fn(mons.len(), |i| {
    coefs[i].to_str() + " " + mons[i].to_str()
  });
  term_strs.connect(" + ")
}

impl<M:Monomial> ToStr
             for PolyOwning<M> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

impl<'self,M:Monomial> ToStr
                   for PolyBorrowing<'self,M> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

impl<'self,M:Monomial> ToStr
                   for PolyBorrowingMons<'self,M> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

pub trait Scalable {

  fn scaled(&self, r: R) -> Self;

  fn scale(&mut self, r: R) -> ();

}

impl<'self,M:Monomial> Scalable
                   for PolyBorrowingMons<'self,M> {
  #[inline]
  fn scaled(&self, r: R) -> PolyBorrowingMons<'self, M> {
    let scaled_coefs = vec::from_fn(self.coefs.len(), |i| r * self.coefs[i]);
    PolyBorrowingMons { coefs: scaled_coefs, mons: self.mons }
  }
  
  #[inline]
  fn scale(&mut self, r: R) -> () {
    for i in range(0, self.coefs.len()) {
      self.coefs[i] *= r;
    }
  }
}

impl<M:Monomial> Scalable
             for PolyOwning<M> {
  #[inline]
  fn scaled(&self, r: R) -> PolyOwning<M> {
    let scaled_coefs = vec::from_fn(self.coefs.len(), |i| r * self.coefs[i]);
    PolyOwning { coefs: scaled_coefs, mons: self.mons.clone() }
  }
  
  #[inline]
  fn scale(&mut self, r: R) -> () {
    for i in range(0, self.coefs.len()) {
      self.coefs[i] *= r;
    }
  }
}


pub fn mul<M:Monomial,P1:Polynomial<M>,P2:Polynomial<M>>(p1: &P1, p2: &P2) -> PolyOwning<M> {
  let n = p1.num_terms() * p2.num_terms();
  let mut mons = vec::with_capacity(n);
  let mut coefs = vec::with_capacity(n);
  p1.each_term(|(c1,m1)| 
    p2.each_term(|(c2,m2)| {
      mons.push(m1 * m2);
      coefs.push(c1 * c2);
    })
  );
  PolyOwning { coefs: coefs, mons: mons }
}


// Convenience function to create polynomials in testing code. Performance critical code
// should use the ::new function implementations instead which will be more efficient.
pub fn poly<M:Monomial>(terms: &[(R,M)]) -> PolyOwning<M> {
  let mut coefs = vec::with_capacity(terms.len());
  let mut mons = vec::with_capacity(terms.len());
  for &(c, ref mon) in terms.iter() {
    coefs.push(c);
    mons.push(mon.clone());
  }
  PolyOwning { coefs: coefs, mons: mons } 
}

pub fn approx_equiv<M:Monomial,P:Polynomial<M>>(p1: &P, p2: &P, tol: R) -> bool {
  let diff = PolyOwning::from_polys_lcomb([(1.,p1),(-1.,p2)]);
  diff.coefs.iter().all(|&c| abs(c) <= tol)
}


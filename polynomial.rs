use extra::treemap::TreeMap;
use std::vec;
use std::num::abs;

use common::*;
use monomial::*;


#[deriving(Eq, Clone)]
pub struct PolyOwning<Mon> {
  coefs: ~[R],
  mons:  ~[Mon]
}

#[deriving(Eq, Clone)]
pub struct PolyBorrowing<'self,Mon> {
  coefs: &'self [R],
  mons: &'self [Mon]
}

#[deriving(Eq, Clone)]
pub struct PolyBorrowingMons<'self,Mon> {
  coefs: ~[R],
  mons: &'self [Mon]
}


impl<Mon:Monomial> PolyOwning<Mon> {
  
  pub fn new(coefs: ~[R], mons: ~[Mon]) -> PolyOwning<Mon> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyOwning { coefs: coefs, mons: mons }
    }
  }
  
  pub fn zero() -> PolyOwning<Mon> {
    let one_mon: Mon = Monomial::one();
    PolyOwning { coefs: ~[0 as R], mons: ~[one_mon] }
  }
  
  pub fn zero_with_capacity(n: uint) -> PolyOwning<Mon> {
    let one_mon: Mon = Monomial::one();
    let mut coefs = vec::with_capacity(n);
    let mut mons = vec::with_capacity(n);
    coefs.push(0 as R);
    mons.push(one_mon);
    PolyOwning { coefs: coefs, mons: mons }
  }
  
  pub fn from_polys_lcomb<P:Polynomial<Mon>>(terms: &[(R,&P)]) -> PolyOwning<Mon> {
    let mut coefs_by_mon: TreeMap<Mon,R> = TreeMap::new();
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
    let mut mons:  ~[Mon] = vec::with_capacity(coefs_by_mon.len());
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

impl<'self,Mon:Monomial> PolyBorrowing<'self,Mon> {
  
  pub fn new(coefs: &'self [R], mons: &'self [Mon]) -> PolyBorrowing<'self,Mon> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyBorrowing { coefs: coefs, mons: mons }
    }
  }
}

impl<'self,Mon:Monomial> PolyBorrowingMons<'self,Mon> {
  
  pub fn new(coefs: ~[R], mons: &'self [Mon]) -> PolyBorrowingMons<'self,Mon> {
    match (coefs.len(), mons.len()) {
      (0,0) => fail!("Empty arrays passed to polynomial constructor."),
      (x,y) if x != y => fail!("Arrays of different lengths passed to polynomial constructor."),
      _ => PolyBorrowingMons { coefs: coefs, mons: mons }
    }
  }
}

pub trait Polynomial<Mon>: ToStr {

  fn num_terms(&self) -> uint;
  
  fn term(&self, n: uint) -> (R,Mon);

  fn foldl_terms<A>(&self, z: A, f: |a: A, term: (R,Mon)| -> A) -> A;
  
  fn foldl_numbered_terms<A>(&self, z: A, f: |a: A, term: (uint,R,Mon)| -> A) -> A;

  fn value_at(&self, x: &[R]) -> R;
  
  fn each_term(&self, f: |term: (R,Mon)| -> ()) ->  ();

  fn canonical_form(&self) -> PolyOwning<Mon>;
  
  fn equiv<P: Polynomial<Mon>>(&self, other: &P) -> bool;

}


fn canonical_form_impl<Mon: Monomial>(coefs: &[R], mons: &[Mon]) -> (~[R],~[Mon]) {
  let mut coefs_by_mon: TreeMap<Mon,R> = TreeMap::new();
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
  let mut mons:  ~[Mon] = vec::with_capacity(coefs_by_mon.len());
  let mut coefs: ~[R] = vec::with_capacity(coefs_by_mon.len());
  for (mon, coef) in coefs_by_mon.iter() {
    if *coef != 0 as R {
      mons.push(mon.clone());
      coefs.push(*coef);
    }
  }
  (coefs, mons)
}


impl<Mon:Monomial> Polynomial<Mon>
               for PolyOwning<Mon> {

  #[inline]
  fn num_terms(&self) -> uint {
    self.mons.len()
  }
  
  #[inline]
  fn term(&self, n: uint) -> (R,Mon) {
    let coef = self.coefs[n];
    let mon = unsafe { self.mons.unsafe_get(n) }; // bounds-checked already in coef access
    (coef,mon)
  }
  
  fn foldl_terms<A>(&self, z: A, f: |a: A, term: (R,Mon)| -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.mons.len()) {
      let term = unsafe { (self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }
  
  fn foldl_numbered_terms<A>(&self, z: A, f: |a: A, term: (uint,R,Mon)| -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.num_terms()) {
      let term = unsafe { (n, self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }
 
  #[inline]
  fn value_at(&self, x: &[R]) -> R {
    self.foldl_terms(0 as R, |sum, (c,m)| sum + c * m.value_at(x))
  }

  #[inline]
  fn each_term(&self, f: |term: (R,Mon)| -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyOwning<Mon> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let m: PolyOwning<Mon> = PolyOwning::zero(); m}
      _ => PolyOwning::new(can_coefs, can_mons)
    }
  }
 
  fn equiv<P: Polynomial<Mon>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }

}

impl<'self,Mon:Monomial> Polynomial<Mon>
                     for PolyBorrowing<'self,Mon> {

  #[inline]
  fn num_terms(&self) -> uint {
    self.mons.len()
  }

  #[inline]
  fn term(&self, n: uint) -> (R,Mon) {
    let coef = self.coefs[n];
    let mon = unsafe { self.mons.unsafe_get(n) }; // bounds-checked already in coef access
    (coef,mon)
  }
  
  fn foldl_terms<A>(&self, z: A, f: |a: A, term: (R,Mon)| -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.mons.len()) {
      let term = unsafe { (self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }
  
  fn foldl_numbered_terms<A>(&self, z: A, f: |a: A, term: (uint,R,Mon)| -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.num_terms()) {
      let term = unsafe { (n, self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }
  
  #[inline]
  fn value_at(&self, x: &[R]) -> R {
    self.foldl_terms(0 as R, |sum, (c,m)| sum + c * m.value_at(x))
  }

  #[inline]
  fn each_term(&self, f: |term: (R,Mon)| -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyOwning<Mon> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let zero: PolyOwning<Mon> = PolyOwning::zero(); zero}
      _ => PolyOwning::new(can_coefs, can_mons)
    }
  }
  
  fn equiv<P: Polynomial<Mon>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }
 
}

impl<'self,Mon:Monomial> Polynomial<Mon>
                     for PolyBorrowingMons<'self,Mon> {

  #[inline]
  fn num_terms(&self) -> uint {
    self.mons.len()
  }

  #[inline]
  fn term(&self, n: uint) -> (R,Mon) {
    let coef = self.coefs[n];
    let mon = unsafe { self.mons.unsafe_get(n) }; // bounds-checked already in coef access
    (coef,mon)
  }
  
  fn foldl_terms<A>(&self, z: A, f: |a: A, term: (R,Mon)| -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.mons.len()) {
      let term = unsafe { (self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }
  
  fn foldl_numbered_terms<A>(&self, z: A, f: |a: A, term: (uint,R,Mon)| -> A) -> A {
    let mut acc_val = z;
    for n in range(0, self.num_terms()) {
      let term = unsafe { (n, self.coefs.unsafe_get(n), self.mons.unsafe_get(n)) };
      acc_val = f(acc_val, term);
    }
    acc_val
  }
  
  #[inline]
  fn value_at(&self, x: &[R]) -> R {
    self.foldl_terms(0 as R, |sum, (c,m)| sum + c * m.value_at(x))
  }

  #[inline]
  fn each_term(&self, f: |term: (R,Mon)| -> ()) ->  () {
    for n in range(0, self.mons.len()) {
      unsafe { f((self.coefs.unsafe_get(n), self.mons.unsafe_get(n))) };
    }
  }

  fn canonical_form(&self) -> PolyOwning<Mon> {
    let (can_coefs, can_mons) = canonical_form_impl(self.coefs, self.mons);
    match can_mons.len() {
      0 => { let zero: PolyOwning<Mon> = PolyOwning::zero(); zero}
      _ => PolyOwning::new(can_coefs, can_mons)
    }
  }
  
  fn equiv<P: Polynomial<Mon>>(&self, other: &P) -> bool {
    let (self_can_coefs, self_can_mons) = canonical_form_impl(self.coefs, self.mons);
    let other_can = other.canonical_form();
    self_can_coefs == other_can.coefs && self_can_mons == other_can.mons 
  }
 
}

// Implement ToStr trait.

fn to_str_impl<Mon: Monomial>(coefs: &[R], mons: &[Mon]) -> ~str {
  let term_strs = vec::from_fn(mons.len(), |i| {
    coefs[i].to_str() + " " + mons[i].to_str()
  });
  term_strs.connect(" + ")
}

impl<Mon:Monomial> ToStr
             for PolyOwning<Mon> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

impl<'self,Mon:Monomial> ToStr
                   for PolyBorrowing<'self,Mon> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

impl<'self,Mon:Monomial> ToStr
                   for PolyBorrowingMons<'self,Mon> {
  fn to_str(&self) -> ~str {
    to_str_impl(self.coefs, self.mons)
  }
}

pub trait Scalable {

  fn scaled(&self, r: R) -> Self;

  fn scale(&mut self, r: R) -> ();

}

impl<'self,Mon:Monomial> Scalable
                   for PolyBorrowingMons<'self,Mon> {
  #[inline]
  fn scaled(&self, r: R) -> PolyBorrowingMons<'self, Mon> {
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

impl<Mon:Monomial> Scalable
             for PolyOwning<Mon> {
  #[inline]
  fn scaled(&self, r: R) -> PolyOwning<Mon> {
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


pub fn mul<Mon:Monomial,P1:Polynomial<Mon>,P2:Polynomial<Mon>>(p1: &P1, p2: &P2) -> PolyOwning<Mon> {
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
pub fn poly<Mon:Monomial>(terms: &[(R,Mon)]) -> PolyOwning<Mon> {
  let mut coefs = vec::with_capacity(terms.len());
  let mut mons = vec::with_capacity(terms.len());
  for &(c, ref mon) in terms.iter() {
    coefs.push(c);
    mons.push(mon.clone());
  }
  PolyOwning { coefs: coefs, mons: mons } 
}

pub fn approx_equiv<Mon:Monomial,P:Polynomial<Mon>>(p1: &P, p2: &P, tol: R) -> bool {
  let diff = PolyOwning::from_polys_lcomb([(1.,p1),(-1.,p2)]);
  diff.coefs.iter().all(|&c| abs(c) <= tol)
}

pub fn approx_equiv_v<Mon:Monomial,P:Polynomial<Mon>>(p1s: &[P], p2s: &[P], tol: R) -> bool {
  p1s.len() == p2s.len() && range(0, p1s.len()).all(|i| approx_equiv(&p1s[i], &p2s[i], tol))
}


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




// Tests


#[test]
#[should_fail]
fn test_1d_owned_constr_mismatched() {
  let one_mon: Mon1d = Monomial::one();
  let bad = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon]);
  println(bad.to_str());
}

#[test]
#[should_fail]
fn test_1d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon1d] = ~[];
  let bad = PolyWithOwnedMons::new(no_Rs, no_mons);
  println(bad.to_str());
}


#[test]
fn test_1d_owned_equiv() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };

  assert!(PolyWithOwnedMons::new(~[1.,2.,1.,3.], ~[one_mon, x_mon, one_mon, x_mon])
  .equiv(&PolyWithOwnedMons::new(~[2.,5.], ~[one_mon, x_mon])))

  let mons = ~[one_mon, x_mon];
  assert!(PolyWithOwnedMons::new(~[1.,2.,1.,3.], ~[one_mon, x_mon, one_mon, x_mon])
  .equiv(&PolyWithBorrowedMons::new(~[2.,5.], mons)))
}  

#[test]
fn test_1d_owned_scaling() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  assert_eq!(PolyWithOwnedMons::new(~[1.,2.], ~[one_mon, x_mon]).scaled(2.),
             PolyWithOwnedMons::new(~[2.,4.], ~[one_mon, x_mon]))
  // scaling by mutation
  let mut p = PolyWithOwnedMons::new(~[1.,2.], ~[one_mon, x_mon]);
  p.scale(3.);
  assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,6.], ~[one_mon, x_mon])));
}  

#[test]
fn test_1d_owned_addition() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let x = PolyWithOwnedMons::new(~[1.], ~[x_mon]);

  assert_eq!(one + x, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,x_mon]));
  assert_eq!(one + one + x, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, one_mon, x_mon]));
  assert_eq!(one + x + x, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, x_mon, x_mon]));
}

#[test]
fn test_1d_owned_multiplication() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let one_plus_x = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, x_mon]);
  let two_plus_x = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, x_mon]);

  assert_eq!(one_plus_x * one_plus_x, 
             PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[one_mon, x_mon, x_mon, x_mon*x_mon]));
  assert_eq!(one_plus_x * two_plus_x,
             PolyWithOwnedMons::new(~[2.,1.,2.,1.], ~[one_mon, x_mon, x_mon, x_mon*x_mon]));
}


#[test]
fn test_1d_owned_canonform() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[x_mon*x_mon, one_mon, x_mon, x_mon]).canonical_form(),
             PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon, x_mon, x_mon*x_mon]));

}

#[test]
fn test_1d_owned_tostr() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.,3.,4.], ~[x_mon*x_mon, one_mon, x_mon, x_mon]).to_str(),
             ~"1 x^2 + 2 x^0 + 3 x^1 + 4 x^1");
}



#[test]
#[should_fail]
fn test_1d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon1d] = ~[Monomial::one()];
  let bad = PolyWithBorrowedMons::new(~[1.,1.], one_mon_singleton);
  println(bad.to_str());
}

#[test]
#[should_fail]
fn test_1d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon1d] = ~[];
  let bad = PolyWithBorrowedMons::new(no_Rs, no_mons);
  println(bad.to_str());
}


#[test]
fn test_1d_borrowed_equiv() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let mons = ~[one_mon, x_mon];
  let one_x_one_x = ~[one_mon, x_mon, one_mon, x_mon];

  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithBorrowedMons::new(~[1.,2.,1.,3.], one_x_one_x)));
  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.,3.], one_x_one_x)));
}  


#[test]
fn test_1d_borrowed_scaling() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let mons = ~[one_mon, x_mon];
  assert_eq!(PolyWithBorrowedMons::new(~[1.,2.], mons).scaled(2.),
             PolyWithBorrowedMons::new(~[2.,4.], mons));
  // scaling by mutation
  let mut p = PolyWithBorrowedMons::new(~[1.,2.], mons);
  p.scale(3.);
  assert!(p.equiv(&PolyWithBorrowedMons::new(~[3.,6.], mons)));
}  


#[test]
fn test_1d_borrowed_addition() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let x_mon_singleton = ~[x_mon];
  let one = PolyWithBorrowedMons::new(~[1.], one_mon_singleton);
  let x = PolyWithBorrowedMons::new(~[1.], x_mon_singleton);

  assert_eq!(one + x, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, x_mon]));
}


/*

#[test]
fn test_2d() {
  // construction
  let one_mon: Mon2d = Monomial::one();
  let x_mon = Mon2d { exps: [Deg(1),Deg(0)] };
  let y_mon = Mon2d { exps: [Deg(0),Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let two = one.scaled(2.0);
  let x = PolyWithOwnedMons::new(~[1.], ~[x_mon]);
  let y = PolyWithOwnedMons::new(~[1.], ~[y_mon]);
  let one_plus_x = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,x_mon]);
  let two_plus_x = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon,x_mon]);
  assert!(one_plus_x.equiv(&one_plus_x));
  assert!(!one_plus_x.equiv(&two_plus_x));

  // scaling
  assert!(one_plus_x.scaled(2.).equiv(&PolyWithOwnedMons::new(~[2.,2.], ~[one_mon, x_mon])));
  assert!(!one_plus_x.scaled(2.).equiv(&PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, x_mon])));
  {
    let mut p = one_plus_x.clone();
    p.scale(3.);
    assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,3.], ~[one_mon, x_mon])));
  }

  // addition and multiplication

  assert!((one + x).equiv(&one_plus_x));
  
  assert!((two + x).equiv(&two_plus_x));

  assert!( ((one + x)*(one + x))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!( ((one + y)*(one + y))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
 
  assert!((one + x.scaled(2.) + x*x)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!((one + y.scaled(2.) + y*y)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );

  assert!( ((one + x)*(one + x)).equiv(&(one + x.scaled(2.) + x*x)) );
  assert!( ((one + y)*(one + y)).equiv(&(one + y.scaled(2.) + y*y)) );

  assert!( ((one + x)*(two + x))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!( ((one + y)*(two + y))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );

  assert!( ((y + x)*(y + x)).equiv(&(y*y + (x*y).scaled(2.) + x*x)) );

  assert_eq!(((one + x)*(one + x)).canonical_form().to_str(),
             ~"1 x^0y^0 + 2 x^1y^0 + 1 x^2y^0");
  assert_eq!(((one + x)*(two + x)).canonical_form().to_str(),
             ~"2 x^0y^0 + 3 x^1y^0 + 1 x^2y^0");
  assert_eq!(((one + y)*(one + y)).canonical_form().to_str(),
             ~"1 x^0y^0 + 2 x^0y^1 + 1 x^0y^2");
  assert_eq!(((one + x)*(two + y)).canonical_form().to_str(),
             ~"2 x^0y^0 + 1 x^0y^1 + 2 x^1y^0 + 1 x^1y^1");
}

#[test]
fn test_3d() {
  // construction
  let one_mon: Mon3d = Monomial::one();
  let x_mon = Mon3d { exps: [Deg(1),Deg(0),Deg(0)] };
  let y_mon = Mon3d { exps: [Deg(0),Deg(1),Deg(0)] };
  let z_mon = Mon3d { exps: [Deg(0),Deg(0),Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let two = one.scaled(2.0);
  let x = PolyWithOwnedMons::new(~[1.], ~[x_mon]);
  let y = PolyWithOwnedMons::new(~[1.], ~[y_mon]);
  let z = PolyWithOwnedMons::new(~[1.], ~[z_mon]);
  let one_plus_x = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,x_mon]);
  let two_plus_x = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon,x_mon]);

  // scaling
  assert!(one_plus_x.scaled(2.).equiv(&PolyWithOwnedMons::new(~[2.,2.], ~[one_mon, x_mon])));
  assert!(!one_plus_x.scaled(2.).equiv(&PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, x_mon])));
  {
    let mut p = one_plus_x.clone();
    p.scale(3.);
    assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,3.], ~[one_mon, x_mon])));
  }

  // addition and multiplication

  assert!((one + x).equiv(&one_plus_x));
  
  assert!((two + x).equiv(&two_plus_x));

  assert!( ((one + x)*(one + x))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!( ((one + y)*(one + y))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
  assert!( ((one + z)*(one + z))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,z_mon,z_mon*z_mon])) );

  assert!( ((one + x)*(one + x)).equiv(&(one + x.scaled(2.) + x*x)) );
  assert!( ((one + y)*(one + y)).equiv(&(one + y.scaled(2.) + y*y)) );
  assert!( ((one + z)*(one + z)).equiv(&(one + z.scaled(2.) + z*z)) );

  assert!((one + x.scaled(2.) + x*x)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!((one + y.scaled(2.) + y*y)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
  assert!((one + z.scaled(2.) + z*z)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,z_mon,z_mon*z_mon])) );

  assert!( ((one + x)*(two + x))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!( ((one + y)*(two + y))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
  assert!( ((one + z)*(two + z))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,z_mon,z_mon*z_mon])) );

  assert!( ((y + x)*(y + x)).equiv(&(y*y + (x*y).scaled(2.) + x*x)) );
  assert!( ((z + x)*(z + x)).equiv(&(z*z + (x*z).scaled(2.) + x*x)) );
  assert!( ((z + y)*(z + y)).equiv(&(z*z + (y*z).scaled(2.) + y*y)) );

  assert_eq!(((one + x)*(one + x)).canonical_form().to_str(),
             ~"1 x^0y^0z^0 + 2 x^1y^0z^0 + 1 x^2y^0z^0");
  assert_eq!(((one + x)*(two + x)).canonical_form().to_str(),
             ~"2 x^0y^0z^0 + 3 x^1y^0z^0 + 1 x^2y^0z^0");
  assert_eq!(((one + y)*(one + y)).canonical_form().to_str(),
             ~"1 x^0y^0z^0 + 2 x^0y^1z^0 + 1 x^0y^2z^0");
  assert_eq!(((one + x)*(two + y)).canonical_form().to_str(),
             ~"2 x^0y^0z^0 + 1 x^0y^1z^0 + 2 x^1y^0z^0 + 1 x^1y^1z^0");
  assert_eq!(((one + y)*(two + z)).canonical_form().to_str(),
             ~"2 x^0y^0z^0 + 1 x^0y^0z^1 + 2 x^0y^1z^0 + 1 x^0y^1z^1");
}

fn some_fn<M: ToStr>(coefs: &[R], mons: &[M]) -> R {
  for i in range(0,mons.len()) {
    println(coefs[i].to_str() + " " + mons[i].to_str());
  }
  1.0
}

#[test]
fn test_4d() {
  // construction
  let one_mon: Mon4d = Monomial::one();
  let x_mon = Mon4d { exps: [Deg(1),Deg(0),Deg(0),Deg(0)] };
  let y_mon = Mon4d { exps: [Deg(0),Deg(1),Deg(0),Deg(0)] };
  let z_mon = Mon4d { exps: [Deg(0),Deg(0),Deg(1),Deg(0)] };
  let t_mon = Mon4d { exps: [Deg(0),Deg(0),Deg(0),Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let two = one.scaled(2.0);
  let x_mon_singleton = ~[x_mon];
  let x = PolyWithOwnedMons::new(~[1.], x_mon_singleton);
  let y = PolyWithOwnedMons::new(~[1.], ~[y_mon]);
  let z = PolyWithOwnedMons::new(~[1.], ~[z_mon]);
  let t = PolyWithOwnedMons::new(~[1.], ~[t_mon]);
  let one_plus_x = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,x_mon]);
  let two_plus_x = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon,x_mon]);


  // TODO: test of borrowing.
  let mons = ~[one_mon, x_mon];
  for i in range(0,2) {
    some_fn(~[i as float, 2.], mons);
  }
  println("Original owned mons: ");
  for mon in mons.iter() {
    println(mon.to_str());
  }
  

  // scaling
  assert!(one_plus_x.scaled(2.).equiv(&PolyWithOwnedMons::new(~[2.,2.], ~[one_mon, x_mon])));
  assert!(!one_plus_x.scaled(2.).equiv(&PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, x_mon])));
  {
    let mut p = one_plus_x.clone();
    p.scale(3.);
    assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,3.], ~[one_mon, x_mon])));
  }

  // addition and multiplication

  assert!((one + x).equiv(&one_plus_x));
  assert!((two + x).equiv(&two_plus_x));

  assert!( ((one + x)*(one + x))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!( ((one + y)*(one + y))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
  assert!( ((one + z)*(one + z))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,z_mon,z_mon*z_mon])) );
  assert!( ((one + t)*(one + t))
             .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,t_mon,t_mon*t_mon])) );

  assert!( ((one + x)*(one + x)).equiv(&(one + x.scaled(2.) + x*x)) );
  assert!( ((one + y)*(one + y)).equiv(&(one + y.scaled(2.) + y*y)) );
  assert!( ((one + z)*(one + z)).equiv(&(one + z.scaled(2.) + z*z)) );
  assert!( ((one + t)*(one + t)).equiv(&(one + t.scaled(2.) + t*t)) );

  assert!((one + x.scaled(2.) + x*x)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!((one + y.scaled(2.) + y*y)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
  assert!((one + z.scaled(2.) + z*z)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,z_mon,z_mon*z_mon])) );
  assert!((one + t.scaled(2.) + t*t)
            .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon,t_mon,t_mon*t_mon])) );

  assert!( ((one + x)*(two + x))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,x_mon,x_mon*x_mon])) );
  assert!( ((one + y)*(two + y))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,y_mon,y_mon*y_mon])) );
  assert!( ((one + z)*(two + z))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,z_mon,z_mon*z_mon])) );
  assert!( ((one + t)*(two + t))
             .equiv(&PolyWithOwnedMons::new(~[2.,3.,1.], ~[one_mon,t_mon,t_mon*t_mon])) );

  assert!( ((y + x)*(y + x)).equiv(&(y*y + (x*y).scaled(2.) + x*x)) );
  assert!( ((z + x)*(z + x)).equiv(&(z*z + (x*z).scaled(2.) + x*x)) );
  assert!( ((z + y)*(z + y)).equiv(&(z*z + (y*z).scaled(2.) + y*y)) );
  assert!( ((z + t)*(z + t)).equiv(&(z*z + (t*z).scaled(2.) + t*t)) );
  assert!( ((t + x)*(t + x)).equiv(&(t*t + (x*t).scaled(2.) + x*x)) );

  assert_eq!(((one + x)*(one + x)).canonical_form().to_str(),
             ~"1 x1^0x2^0x3^0x4^0 + 2 x1^1x2^0x3^0x4^0 + 1 x1^2x2^0x3^0x4^0");
  assert_eq!(((one + x)*(two + x)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 3 x1^1x2^0x3^0x4^0 + 1 x1^2x2^0x3^0x4^0");
  assert_eq!(((one + y)*(one + y)).canonical_form().to_str(),
             ~"1 x1^0x2^0x3^0x4^0 + 2 x1^0x2^1x3^0x4^0 + 1 x1^0x2^2x3^0x4^0");
  assert_eq!(((one + x)*(two + y)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 1 x1^0x2^1x3^0x4^0 + 2 x1^1x2^0x3^0x4^0 + 1 x1^1x2^1x3^0x4^0");
  assert_eq!(((one + y)*(two + z)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 1 x1^0x2^0x3^1x4^0 + 2 x1^0x2^1x3^0x4^0 + 1 x1^0x2^1x3^1x4^0");

  assert_eq!(((one + t)*(two + t)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 3 x1^0x2^0x3^0x4^1 + 1 x1^0x2^0x3^0x4^2");
  assert_eq!(((one + x)*(two + t)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 1 x1^0x2^0x3^0x4^1 + 2 x1^1x2^0x3^0x4^0 + 1 x1^1x2^0x3^0x4^1");
  assert_eq!(((one + y)*(two + t)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 1 x1^0x2^0x3^0x4^1 + 2 x1^0x2^1x3^0x4^0 + 1 x1^0x2^1x3^0x4^1");
  assert_eq!(((one + z)*(two + t)).canonical_form().to_str(),
             ~"2 x1^0x2^0x3^0x4^0 + 1 x1^0x2^0x3^0x4^1 + 2 x1^0x2^0x3^1x4^0 + 1 x1^0x2^0x3^1x4^1");
}
*/

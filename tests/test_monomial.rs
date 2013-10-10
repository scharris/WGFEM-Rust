use common::*;
use monomial;
use monomial::{Monomial, Mon1d, Mon2d, Mon3d, Mon4d, MaxMonDeg};

fn test_domain_dims() {
  assert_eq!(Monomial::domain_space_dims(None::<Mon1d>), 1u);
  assert_eq!(Monomial::domain_space_dims(None::<Mon2d>), 2u);
  assert_eq!(Monomial::domain_space_dims(None::<Mon3d>), 3u);
  assert_eq!(Monomial::domain_space_dims(None::<Mon4d>), 4u);
  
  assert_eq!(monomial::domain_space_dims::<Mon1d>(), 1u);
  assert_eq!(monomial::domain_space_dims::<Mon2d>(), 2u);
  assert_eq!(monomial::domain_space_dims::<Mon3d>(), 3u);
  assert_eq!(monomial::domain_space_dims::<Mon4d>(), 4u);
}

#[test]
fn test_value_at_1d() {
  let one = Mon1d { exps: [Deg(0)] };
  let x = Mon1d { exps: [Deg(1)] };
  let x2 = Mon1d { exps: [Deg(2)] };
  let x3 = Mon1d { exps: [Deg(3)] };
  assert_eq!(one.value_at(&[2.]), 1.);
  assert_eq!(x.value_at(&[2.]), 2.);
  assert_eq!(x2.value_at(&[2.]), 4.);
  assert_eq!(x3.value_at(&[2.]), 8.);
}
#[test]
fn test_value_at_2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  assert_eq!(one.value_at(&[2.,3.]), 1.);
  
  let x1y2 = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(x1y2.value_at(&[2.,3.]), 18.);
  
  let x3y1 = Mon2d { exps: [Deg(3), Deg(1)] };
  assert_eq!(x3y1.value_at(&[2.,3.]), 24.);
}
#[test]
fn test_value_at_3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at(&[2.,3.,4.]), 1.);
  
  let x1y2z3 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(x1y2z3.value_at(&[2.,3.,4.]), 18.*64.);
  
  let x3y1z2 = Mon3d { exps: [Deg(3), Deg(1), Deg(2)] };
  assert_eq!(x3y1z2.value_at(&[2.,3.,4.]), 24.*16.);
}
#[test]
fn test_value_at_4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at(&[2.,3.,4.,5.]), 1.);
  
  let x1y2z3t4 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(x1y2z3t4.value_at(&[2.,3.,4.,5.]), 18.*64.*5.*5.*5.*5.);
  
  let x3y1z2t4 = Mon4d { exps: [Deg(3), Deg(1), Deg(2), Deg(4)] };
  assert_eq!(x3y1z2t4.value_at(&[2.,3.,4.,5.]), 24.*16.*5.*5.*5.*5.);
}

#[test]
fn test_value_at_for_origin_1d() {
  let one = Mon1d { exps: [Deg(0)] };
  let x = Mon1d { exps: [Deg(1)] };
  let x2 = Mon1d { exps: [Deg(2)] };
  let x3 = Mon1d { exps: [Deg(3)] };

  assert_eq!(one.value_at_for_origin(&[3.],&[1.]), 1.);
  assert_eq!(x.value_at_for_origin(&[3.],&[1.]), 2.);
  assert_eq!(x2.value_at_for_origin(&[3.],&[1.]), 4.);
  assert_eq!(x3.value_at_for_origin(&[3.],&[1.]), 8.);
}
#[test]
fn test_value_at_for_origin_2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  assert_eq!(one.value_at_for_origin(&[2.,3.],&[12.,13.]), 1.);
  
  let x1y2 = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(x1y2.value_at_for_origin(&[2.,3.],&[-1.,-2.]), 3. * 25.);
  
  let x3y1 = Mon2d { exps: [Deg(3), Deg(1)] };
  assert_eq!(x3y1.value_at_for_origin(&[2.,3.],&[-1.,1.]), 27. * 2.);
}
#[test]
fn test_value_at_for_origin_3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at_for_origin(&[2.,3.,4.],&[1.,0.,2.]), 1.);
  
  let x1y2z3 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(x1y2z3.value_at_for_origin(&[2.,3.,4.],&[3.,1.,2.]), -1. * 4. * 8.);
  
  let x3y1z2 = Mon3d { exps: [Deg(3), Deg(1), Deg(2)] };
  assert_eq!(x3y1z2.value_at_for_origin(&[2.,3.,4.],&[3.,1.,2.]), -1. * 2. * 4.);
}
#[test]
fn test_value_at_for_origin_4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  assert_eq!(one.value_at_for_origin(&[2.,3.,4.,5.],&[1.,20.,30.,40.]), 1.);
  
  let x1y2z3t4 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(x1y2z3t4.value_at_for_origin(&[2.,3.,4.,5.],&[0.,5.,1.,2.]), 2. * 4. * 27. * 3.*3.*3.*3.);
  
  let x3y1z2t4 = Mon4d { exps: [Deg(3), Deg(1), Deg(2), Deg(4)] };
  assert_eq!(x3y1z2t4.value_at_for_origin(&[2.,3.,4.,5.],&[0.,5.,1.,2.]), 8. * (-2.) * 9. * 3.*3.*3.*3.);
}


#[test]
fn test_exp_1d() {
  let m = Mon1d { exps: [Deg(1)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
}

#[test]
#[should_fail]
fn test_exp_fail_1d() {
  let one: Mon1d = Monomial::one();
  one.exp(Dim(1)); 
}

#[test]
fn test_exp_2d() {
  let m = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
  assert_eq!(m.exp(Dim(1)), Deg(2));
}
#[test]
#[should_fail]
fn test_exp_fail_2d() {
  let one: Mon2d = Monomial::one();
  one.exp(Dim(2)); 
}
#[test]
fn test_exp_3d() {
  let m = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
  assert_eq!(m.exp(Dim(1)), Deg(2));
  assert_eq!(m.exp(Dim(2)), Deg(3));
}
#[test]
#[should_fail]
fn test_exp_fail_3d() {
  let one: Mon3d = Monomial::one();
  one.exp(Dim(3)); 
}
#[test]
fn test_exp_4d() {
  let m = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(m.exp(Dim(0)), Deg(1));
  assert_eq!(m.exp(Dim(1)), Deg(2));
  assert_eq!(m.exp(Dim(2)), Deg(3));
  assert_eq!(m.exp(Dim(3)), Deg(4));
}
#[test]
#[should_fail]
fn test_exp_fail_4d() {
  let one: Mon4d = Monomial::one();
  one.exp(Dim(4)); 
}

#[test]
fn test_equality_1d() {
  let m1 = Mon1d { exps: [Deg(1)] };
  let m2 = Mon1d { exps: [Deg(1)] };
  assert_eq!(m1, m2);

  let m3 = Mon1d { exps: [Deg(0)] };
  assert!(m1 != m3);
  assert!(!(m1 == m3));
}
#[test]
fn test_equality_2d() {
  let m1 = Mon2d { exps: [Deg(1), Deg(2)] };
  let m2 = Mon2d { exps: [Deg(1), Deg(2)] };
  assert_eq!(m1, m2);

  let m3 = Mon2d { exps: [Deg(1), Deg(0)] };
  let m4 = Mon2d { exps: [Deg(0), Deg(2)] };
  assert!(m1 != m3);
  assert!(m1 != m4);
  assert!(!(m1 == m3));
}
#[test]
fn test_equality_3d() {
  let m1 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  let m2 = Mon3d { exps: [Deg(1), Deg(2), Deg(3)] };
  assert_eq!(m1, m2);

  let m3 = Mon3d { exps: [Deg(0), Deg(2), Deg(3)] };
  let m4 = Mon3d { exps: [Deg(1), Deg(0), Deg(3)] };
  let m5 = Mon3d { exps: [Deg(1), Deg(2), Deg(0)] };
  assert!(m1 != m3);
  assert!(m1 != m4);
  assert!(m1 != m5);
  assert!(!(m1 == m3));
}
#[test]
fn test_equality_4d() {
  let m1 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  let m2 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(4)] };
  assert_eq!(m1, m2);

  let m3 = Mon4d { exps: [Deg(0), Deg(2), Deg(3), Deg(4)] };
  let m4 = Mon4d { exps: [Deg(1), Deg(0), Deg(3), Deg(4)] };
  let m5 = Mon4d { exps: [Deg(1), Deg(2), Deg(0), Deg(4)] };
  let m6 = Mon4d { exps: [Deg(1), Deg(2), Deg(3), Deg(0)] };
  assert!(m1 != m3);
  assert!(m1 != m4);
  assert!(m1 != m5);
  assert!(m1 != m6);
  assert!(!(m1 == m3));
}

#[test]
fn test_ord_1d() {
  let x = Mon1d { exps: [Deg(1)] };
  let one: Mon1d = Monomial::one();
  assert!(one < x);
  assert!(!(x < one));
  assert!(x > one);

  let x2 = Mon1d { exps: [Deg(2)] };
  assert!(x2 > x);
  assert!(x < x2);
  assert!(!(x > x2));
  assert!(!(x2 < x));
}
#[test]
fn test_ord_2d() {
  let one: Mon2d = Monomial::one();
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  assert!(one < y);
  assert!(x > one);
  assert!(!(x < one));

  assert!(y < x);
  assert!(y < x*y);
  assert!(x < x*y);
  assert!(x*y > x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  assert!(x*x > x);
  assert!(x < x*x);
  assert!(!(x > x*x));
  assert!(!(x*x < x));

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);
}
#[test]
fn test_ord_3d() {
  let one: Mon3d = Monomial::one();
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  assert!(one < y);
  assert!(x > one);
  assert!(!(x < one));

  assert!(y < x);
  assert!(y < x*y);
  assert!(x < x*y);
  assert!(x*y > x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  assert!(x*x > x);
  assert!(x < x*x);
  assert!(!(x > x*x));
  assert!(!(x*x < x));

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);

  assert!(z*z < x*y);
  assert!(y*z < x);
  assert!(x*y*z > x);
  assert!(z <= z);
}
#[test]
fn test_ord_4d() {
  let one: Mon4d = Monomial::one();
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  assert!(one < y);
  assert!(x > one);
  assert!(!(x < one));

  assert!(y < x);
  assert!(y < x*y);
  assert!(x < x*y);
  assert!(x*y > x);
  assert!(!(x < y));
  assert!(y <= x);
  assert!(x <= x);

  assert!(one.cmp(&x) == Less);
  assert!(x.cmp(&one) == Greater);
  assert!(x.cmp(&x) == Equal);

  assert!(x*x > x);
  assert!(x < x*x);
  assert!(!(x > x*x));
  assert!(!(x*x < x));

  assert!(y.cmp(&x) == Less);
  assert!(x.cmp(&y) == Greater);
  assert!(x.cmp(&x) == Equal);
  assert!(y.cmp(&y) == Equal);

  assert!(z*z < x*y);
  assert!(y*z < x);
  assert!(x*y*z > x);
  assert!(z <= z);

  assert!(y*z*t*t < x*y*z);
  assert!(t < z);
  assert!(t < y);
  assert!(y*z*t > t);
  assert!(t <= t);
}

#[test]
fn test_clone_1d() {
  let x = Mon1d { exps: [Deg(1)] };
  assert_eq!(x, x.clone());
  let one = Mon1d { exps: [Deg(0)] };
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}
#[test]
fn test_clone_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  assert_eq!(x, x.clone());
  let one: Mon2d = Monomial::one();
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}
#[test]
fn test_clone_3d() {
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  assert_eq!(x, x.clone());
  let one: Mon3d = Monomial::one();
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}
#[test]
fn test_clone_4d() {
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  assert_eq!(x, x.clone());
  let one: Mon4d = Monomial::one();
  assert_eq!(one, one.clone());
  assert!(!(x.clone() == one));
}


#[test]
fn test_hash_1d() {
  use std::hashmap::HashSet;
  let x = Mon1d { exps: [Deg(1)] };
  let one: Mon1d = Monomial::one();
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  m.insert(one);
  assert!(m.contains(&one));
}
#[test]
fn test_hash_2d() {
  use std::hashmap::HashSet;
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one: Mon2d = Monomial::one();
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  assert!(!m.contains(&y));
  m.insert(one);
  assert!(m.contains(&one));
  m.insert(y);
  assert!(m.contains(&y));
}
#[test]
fn test_hash_3d() {
  use std::hashmap::HashSet;
  let one: Mon3d = Monomial::one();
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  assert!(!m.contains(&y));
  assert!(!m.contains(&z));
  m.insert(one);
  assert!(m.contains(&one));
  m.insert(y);
  assert!(m.contains(&y));
  m.insert(z);
  assert!(m.contains(&z));
}
#[test]
fn test_hash_4d() {
  use std::hashmap::HashSet;
  let one: Mon4d = Monomial::one();
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let mut m = HashSet::new();
  m.insert(x);
  assert!(m.contains(&x));
  assert!(!m.contains(&one));
  assert!(!m.contains(&y));
  assert!(!m.contains(&z));
  assert!(!m.contains(&t));
  m.insert(one);
  assert!(m.contains(&one));
  m.insert(y);
  assert!(m.contains(&y));
  m.insert(z);
  assert!(m.contains(&z));
  m.insert(t);
  assert!(m.contains(&t));
}


#[test]
fn test_mul_1d() {
  let one: Mon1d = Monomial::one();
  let x = Mon1d { exps: [Deg(1)] };
  assert_eq!(x * x, Mon1d { exps: [Deg(2)] });
  assert_eq!(x * one, x);
}
fn test_mul_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one: Mon2d = Monomial::one();
  assert_eq!(x * y, Mon2d { exps: [Deg(1), Deg(1)] });
  assert_eq!(x * one, x);
  assert_eq!(one * y, y);
}

#[test]
fn test_mul_3d() {
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one: Mon3d = Monomial::one();
  assert_eq!(x * y, Mon3d { exps: [Deg(1), Deg(1), Deg(0)] });
  assert_eq!(x * one, x);
  assert_eq!(one * y, y);
  assert_eq!(x*z, Mon3d { exps: [Deg(1), Deg(0), Deg(1)] });
  assert_eq!(x*y, Mon3d { exps: [Deg(1), Deg(1), Deg(0)] });
  assert_eq!(y*z, Mon3d { exps: [Deg(0), Deg(1), Deg(1)] });
}
#[test]
fn test_mul_4d() {
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one: Mon4d = Monomial::one();
  assert_eq!(x * y * z * t, Mon4d { exps: [Deg(1), Deg(1), Deg(1), Deg(1)] });
  assert_eq!(x * one, x);
  assert_eq!(x * one, x);
  assert_eq!(z * one, z);
  assert_eq!(t * one, t);
  assert_eq!(x*z, Mon4d { exps: [Deg(1), Deg(0), Deg(1), Deg(0)] });
  assert_eq!(x*y, Mon4d { exps: [Deg(1), Deg(1), Deg(0), Deg(0)] });
  assert_eq!(y*z, Mon4d { exps: [Deg(0), Deg(1), Deg(1), Deg(0)] });
  assert_eq!(t*z, Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(1)] });
  assert_eq!(t*y, Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(1)] });
  assert_eq!(t*x, Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(1)] });
}


#[test]
fn test_map_exp_1d() {
  let x = Mon1d { exps: [Deg(1)] };
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
}
#[test]
fn test_map_exp_2d() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), x*x*y);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
}
#[test]
fn test_map_exp_3d() {
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(x.map_exp(Dim(2), |e| Deg(*e+1)), x*z);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), x*x*y);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
  assert_eq!(y.map_exp(Dim(2), |e| Deg(*e+2)), y*z*z);
  assert_eq!(z.map_exp(Dim(0), |e| Deg(*e+3)), x*x*x*z);
  assert_eq!(z.map_exp(Dim(1), |e| Deg(*e+3)), y*y*y*z);
  assert_eq!(z.map_exp(Dim(2), |e| Deg(*e+3)), z*z*z*z);
}
#[test]
fn test_map_exp_4d() {
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  assert_eq!(x.map_exp(Dim(0), |e| Deg(*e+1)), x*x);
  assert_eq!(x.map_exp(Dim(1), |e| Deg(*e+1)), x*y);
  assert_eq!(x.map_exp(Dim(2), |e| Deg(*e+1)), x*z);
  assert_eq!(x.map_exp(Dim(3), |e| Deg(*e+1)), x*t);
  assert_eq!(y.map_exp(Dim(0), |e| Deg(*e+2)), y*x*x);
  assert_eq!(y.map_exp(Dim(1), |e| Deg(*e+2)), y*y*y);
  assert_eq!(y.map_exp(Dim(2), |e| Deg(*e+2)), y*z*z);
  assert_eq!(y.map_exp(Dim(3), |e| Deg(*e+2)), y*t*t);
  assert_eq!(z.map_exp(Dim(0), |e| Deg(*e+3)), z*x*x*x);
  assert_eq!(z.map_exp(Dim(1), |e| Deg(*e+3)), z*y*y*y);
  assert_eq!(z.map_exp(Dim(2), |e| Deg(*e+3)), z*z*z*z);
  assert_eq!(z.map_exp(Dim(3), |e| Deg(*e+3)), z*t*t*t);
  assert_eq!(t.map_exp(Dim(0), |e| Deg(*e+1)), t*x);
  assert_eq!(t.map_exp(Dim(1), |e| Deg(*e+1)), t*y);
  assert_eq!(t.map_exp(Dim(2), |e| Deg(*e+1)), t*z);
  assert_eq!(t.map_exp(Dim(3), |e| Deg(*e+1)), t*t);
}

#[test]
fn test_num_mons_with_deg_lim() {
  assert_eq!(monomial::num_mons_with_deg_lim(MaxMonDeg(2), 2), 6);

  let mons_4d_deg_le_5: ~[Mon4d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(5));
  assert_eq!(mons_4d_deg_le_5.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(5), 4));

  assert_eq!(monomial::num_mons_with_deg_lim(MaxMonDeg(1), 3), 4);

  assert_eq!(monomial::num_mons_with_deg_lim(MaxMonDeg(0), 2), 1);
}

#[test]
fn test_mons_with_deg_lim_asc_2d() {
  let one: Mon2d = Monomial::one();
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  
  let mons_2d_deg_le_0: ~[Mon2d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(0));
  assert_eq!(&mons_2d_deg_le_0, &~[one]);
  assert_eq!(mons_2d_deg_le_0.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(0), 2));
  
  let mons_2d_deg_le_1: ~[Mon2d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(1));
  assert_eq!(&mons_2d_deg_le_1, &~[one, y, x]);
  assert_eq!(mons_2d_deg_le_1.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(1), 2));
 
  let mons_2d_deg_le_2: ~[Mon2d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(2));
  assert_eq!(&mons_2d_deg_le_2, &~[one, y, y*y, x, x*y, x*x]);
  assert_eq!(mons_2d_deg_le_2.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(2), 2));
}

#[test]
fn test_mons_with_deg_lim_asc_3d() {
  let one: Mon3d = Monomial::one();
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  
  let mons_3d_deg_le_0: ~[Mon3d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(0));
  assert_eq!(&mons_3d_deg_le_0, &~[one]);
  assert_eq!(mons_3d_deg_le_0.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(0), 3));
  
  let mons_3d_deg_le_1: ~[Mon3d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(1));
  assert_eq!(&mons_3d_deg_le_1, &~[one, z, y, x]);
  assert_eq!(mons_3d_deg_le_1.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(1), 3));
 
  let mons_3d_deg_le_2: ~[Mon3d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(2));
  assert_eq!(&mons_3d_deg_le_2, &~[one, z, z*z, y, y*z, y*y, x, x*z, x*y, x*x]);
  assert_eq!(mons_3d_deg_le_2.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(2), 3));
}


#[test]
fn test_mons_with_deg_lim_asc_4d() {
  let one: Mon4d = Monomial::one();
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  
  let mons_4d_deg_le_0: ~[Mon4d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(0));
  assert_eq!(&mons_4d_deg_le_0, &~[one]);
  assert_eq!(mons_4d_deg_le_0.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(0), 4));
  
  let mons_4d_deg_le_1: ~[Mon4d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(1));
  assert_eq!(&mons_4d_deg_le_1, &~[one, t, z, y, x]);
  assert_eq!(mons_4d_deg_le_1.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(1), 4));
 
  let mons_4d_deg_le_2: ~[Mon4d] = Monomial::mons_with_deg_lim_asc(MaxMonDeg(2));
  assert_eq!(&mons_4d_deg_le_2, 
             &~[one, t, t*t, z, z*t, z*z, y, y*t, y*z, y*y, x, x*t, x*z, x*y, x*x]);
  assert_eq!(mons_4d_deg_le_2.len(), monomial::num_mons_with_deg_lim(MaxMonDeg(2), 4));
}


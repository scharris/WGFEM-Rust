use common::*;
use monomial::{Monomial, Mon1d, Mon2d, Mon3d, Mon4d};
use polynomial::*;

#[test]
#[should_fail]
fn test_1d_owned_constr_mismatched() {
  let one_mon: Mon1d = Monomial::one();
  PolyOwning::new(~[1.,1.], ~[one_mon]);
}

#[test]
#[should_fail]
fn test_1d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon1d] = ~[];
  PolyOwning::new(no_Rs, no_mons);
}


#[test]
fn test_1d_owned_equiv() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };

  assert!(PolyOwning::new(~[1.,2.,1.,3.], ~[one_mon, x_mon, one_mon, x_mon])
  .equiv(&PolyOwning::new(~[2.,5.], ~[one_mon, x_mon])))

  let mons = ~[one_mon, x_mon];
  assert!(PolyOwning::new(~[1.,2.,1.,3.], ~[one_mon, x_mon, one_mon, x_mon])
  .equiv(&PolyBorrowingMons::new(~[2.,5.], mons)))
}  

fn test_1d_owned_scaling() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  assert_eq!(PolyOwning::new(~[1.,2.], ~[one_mon, x_mon]).scaled(2.),
             PolyOwning::new(~[2.,4.], ~[one_mon, x_mon]))
}

#[test]
fn test_1d_owned_lcomb() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let one = PolyOwning::new(~[1.], ~[one_mon]);
  let x = PolyOwning::new(~[1.], ~[x_mon]);
  let x2 = PolyOwning::new(~[1.], ~[x_mon*x_mon]);

  assert_eq!(PolyOwning::from_polys_lcomb([(2.,&one), (3.,&x), (2.5,&x), (2.,&x2), (1.,&one)]), poly([(3.,one_mon),(5.5,x_mon),(2.,x_mon*x_mon)]));
}


#[test]
fn test_1d_owned_multiplication() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let one_plus_x = PolyOwning::new(~[1.,1.], ~[one_mon, x_mon]);
  let two_plus_x = PolyOwning::new(~[2.,1.], ~[one_mon, x_mon]);

  assert_eq!(mul_polys(&one_plus_x,  &one_plus_x), 
             PolyOwning::new(~[1.,1.,1.,1.], ~[one_mon, x_mon, x_mon, x_mon*x_mon]));
  assert_eq!(mul_polys(&one_plus_x, &two_plus_x),
             PolyOwning::new(~[2.,1.,2.,1.], ~[one_mon, x_mon, x_mon, x_mon*x_mon]));
}


#[test]
fn test_1d_owned_canonform() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,1.,1.,1.], ~[x_mon*x_mon, one_mon, x_mon, x_mon]).canonical_form(),
             PolyOwning::new(~[1.,2.,1.], ~[one_mon, x_mon, x_mon*x_mon]));

}

#[test]
fn test_1d_owned_tostr() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,2.,3.,4.], ~[x_mon*x_mon, one_mon, x_mon, x_mon]).to_str(),
             ~"1 x^2 + 2 x^0 + 3 x^1 + 4 x^1");
}



#[test]
#[should_fail]
fn test_1d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon1d] = ~[Monomial::one()];
  PolyBorrowingMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_1d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon1d] = ~[];
  PolyBorrowingMons::new(no_Rs, no_mons);
}


#[test]
fn test_1d_borrowed_equiv() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let mons = ~[one_mon, x_mon];
  let one_x_one_x = ~[one_mon, x_mon, one_mon, x_mon];

  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyBorrowingMons::new(~[1.,2.,1.,3.], one_x_one_x)));
  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyOwning::new(~[1.,2.,1.,3.], one_x_one_x)));
}  


#[test]
fn test_1d_borrowed_scaling() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let mons = ~[one_mon, x_mon];
  assert_eq!(PolyBorrowingMons::new(~[1.,2.], mons).scaled(2.),
             PolyBorrowingMons::new(~[2.,4.], mons));
}


#[test]
fn test_1d_borrowed_lcomb() {
  let one_mon: Mon1d = Monomial::one();
  let x_mon = Mon1d { exps: [Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let x_mon_singleton = ~[x_mon];
  let x2_mon_singleton = ~[x_mon*x_mon];
  let one = PolyBorrowingMons::new(~[1.], one_mon_singleton);
  let x = PolyBorrowingMons::new(~[1.], x_mon_singleton);
  let x2 = PolyBorrowingMons::new(~[1.], x2_mon_singleton); 

  assert_eq!(PolyOwning::from_polys_lcomb([(2.,&one), (3.,&x), (2.5,&x), (2.,&x2), (1.,&one)]), poly([(3.,one_mon),(5.5,x_mon),(2.,x_mon*x_mon)]));
}

#[test]
#[should_fail]
fn test_2d_owned_constr_mismatched() {
  let one_mon: Mon2d = Monomial::one();
  PolyOwning::new(~[1.,1.], ~[one_mon]);
}


#[test]
#[should_fail]
fn test_2d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon2d] = ~[];
  PolyOwning::new(no_Rs, no_mons);
}


#[test]
fn test_2d_owned_equiv() {
  let one: Mon2d = Monomial::one();
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert!(poly([(1.,one),(2.,y), (1.,one), (3.,y)])
            .equiv(&poly([(2.,one),(5.,y)])));

  let mons = ~[one, y];
  assert!(poly([(1.,one),(2.,y), (1.,one), (3.,y)])
            .equiv(&PolyBorrowingMons::new(~[2.,5.], mons)));
}  

#[test]
fn test_2d_owned_approx_equiv() {
  let one: Mon2d = Monomial::one();
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert!(approx_equiv(&poly([(1.,one),(2.0000001,y), (1.,one), (3.,y)]), &poly([(2.,one),(5.,y)]), 0.001));

  assert!(!approx_equiv(&poly([(1.,one), (2.01,y), (1.,one), (3.,y)]), &poly([(2.,one), (5.,y)]), 1e-3));
}

#[test]
fn test_2d_owned_scaling() {
  let one: Mon2d = Monomial::one();
  let y = Mon2d { exps: [Deg(1), Deg(0)] };

  assert_eq!(PolyOwning::new(~[1.,2.], ~[one, y]).scaled(2.),
             PolyOwning::new(~[2.,4.], ~[one, y]))
}

#[test]
fn test_2d_owned_lcomb() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one = PolyOwning::new(~[1.], ~[one_mon]);
  let y = PolyOwning::new(~[1.], ~[y_mon]);

  assert_eq!(PolyOwning::from_polys_lcomb([(2.,&one), (3.,&y), (2.5,&y), (1.,&one)]), poly([(3.,one_mon),(5.5,y_mon)]));
}

#[test]
fn test_2d_owned_multiplication() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one_plus_y = PolyOwning::new(~[1.,1.], ~[one_mon, y_mon]);
  let two_plus_y = PolyOwning::new(~[2.,1.], ~[one_mon, y_mon]);

  assert_eq!(mul_polys(&one_plus_y, &one_plus_y), 
             PolyOwning::new(~[1.,1.,1.,1.], ~[one_mon, y_mon, y_mon, y_mon*y_mon]));
  assert_eq!(mul_polys(&one_plus_y, &two_plus_y),
             PolyOwning::new(~[2.,1.,2.,1.], ~[one_mon, y_mon, y_mon, y_mon*y_mon]));
}


#[test]
fn test_2d_owned_canonform() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,1.,1.,1.], ~[y_mon*y_mon, one_mon, y_mon, y_mon]).canonical_form(),
             PolyOwning::new(~[1.,2.,1.], ~[one_mon, y_mon, y_mon*y_mon]));

}

#[test]
fn test_2d_owned_tostr() {
  let x_mon = Mon2d { exps: [Deg(1), Deg(0)] };
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,2.,3.,4.], ~[x_mon*x_mon, x_mon*y_mon, x_mon, x_mon]).to_str(),
             ~"1 x^2y^0 + 2 x^1y^1 + 3 x^1y^0 + 4 x^1y^0");
}



#[test]
#[should_fail]
fn test_2d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon2d] = ~[Monomial::one()];
  PolyBorrowingMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_2d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon2d] = ~[];
  PolyBorrowingMons::new(no_Rs, no_mons);
}


#[test]
fn test_2d_borrowed_equiv() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let mons = ~[one_mon, y_mon];
  let one_y_one_y = ~[one_mon, y_mon, one_mon, y_mon];

  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyBorrowingMons::new(~[1.,2.,1.,3.], one_y_one_y)));
  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyOwning::new(~[1.,2.,1.,3.], one_y_one_y)));
}  


#[test]
fn test_2d_borrowed_scaling() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let mons = ~[one_mon, y_mon];
  assert_eq!(PolyBorrowingMons::new(~[1.,2.], mons).scaled(2.),
             PolyBorrowingMons::new(~[2.,4.], mons));
}


#[test]
fn test_2d_borrowed_lcomb() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let y_mon_singleton = ~[y_mon];
  let one = PolyBorrowingMons::new(~[1.], one_mon_singleton);
  let y = PolyBorrowingMons::new(~[1.], y_mon_singleton);

  assert_eq!(PolyOwning::from_polys_lcomb([(2.,&one), (3.,&y), (2.5,&y), (1.,&one)]), poly([(3.,one_mon),(5.5,y_mon)]));
}

#[test]
fn test_2d_borrowed_multiplication() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one_y_mons = ~[one_mon, y_mon];

  let one_plus_y = PolyBorrowingMons::new(~[1.,1.], one_y_mons);
  let two_plus_y = PolyBorrowingMons::new(~[2.,1.], one_y_mons);

  assert_eq!(mul_polys(&one_plus_y, &one_plus_y), 
             PolyOwning::new(~[1.,1.,1.,1.], ~[one_mon, y_mon, y_mon, y_mon*y_mon]));
  assert_eq!(mul_polys(&one_plus_y, &two_plus_y),
             PolyOwning::new(~[2.,1.,2.,1.], ~[one_mon, y_mon, y_mon, y_mon*y_mon]));
}


#[test]
#[should_fail]
fn test_3d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon3d] = ~[Monomial::one()];
  PolyBorrowingMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_3d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon3d] = ~[];
  PolyBorrowingMons::new(no_Rs, no_mons);
}

#[test]
fn test_3d_borrowed_equiv() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, z_mon];
  let one_z_one_z = ~[one_mon, z_mon, one_mon, z_mon];

  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyBorrowingMons::new(~[1.,2.,1.,3.], one_z_one_z)));
  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyOwning::new(~[1.,2.,1.,3.], one_z_one_z)));
}  


#[test]
fn test_3d_borrowed_scaling() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, z_mon];
  assert_eq!(PolyBorrowingMons::new(~[1.,2.], mons).scaled(2.),
             PolyBorrowingMons::new(~[2.,4.], mons));
}


#[test]
fn test_3d_borrowed_lcomb() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let z_mon_singleton = ~[z_mon];
  let one = PolyBorrowingMons::new(~[1.], one_mon_singleton);
  let z = PolyBorrowingMons::new(~[1.], z_mon_singleton);

  assert_eq!(PolyOwning::from_polys_lcomb([(-2.,&one), (3.,&z), (-2.5,&z), (1.,&one)]), poly([(-1.,one_mon),(0.5,z_mon)]));
}

#[test]
#[should_fail]
fn test_3d_owned_constr_mismatched() {
  let one_mon: Mon3d = Monomial::one();
  PolyOwning::new(~[1.,1.], ~[one_mon]);
}

#[test]
#[should_fail]
fn test_3d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon3d] = ~[];
  PolyOwning::new(no_Rs, no_mons);
}


#[test]
fn test_3d_owned_equiv() {
  let one: Mon3d = Monomial::one();
  let z = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };

  assert!(poly([(1.,one),(2.,z), (1.,one), (3.,z)])
            .equiv(&poly([(2.,one),(5.,z)])));

  let mons = ~[one, z];
  assert!(poly([(1.,one),(2.,z), (1.,one), (3.,z)])
            .equiv(&PolyBorrowingMons::new(~[2.,5.], mons)));
}  

#[test]
fn test_3d_owned_scaling() {
  let one: Mon3d = Monomial::one();
  let z = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };

  assert_eq!(PolyOwning::new(~[1.,2.], ~[one, z]).scaled(2.),
             PolyOwning::new(~[2.,4.], ~[one, z]))
}

#[test]
fn test_3d_owned_lcomb() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one = PolyOwning::new(~[1.], ~[one_mon]);
  let z = PolyOwning::new(~[1.], ~[z_mon]);

  assert_eq!(PolyOwning::from_polys_lcomb([(-2.,&one), (3.,&z), (-2.5,&z), (1.,&one)]), poly([(-1.,one_mon),(0.5,z_mon)]));
}

#[test]
fn test_3d_owned_multiplication() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_plus_z = PolyOwning::new(~[1.,1.], ~[one_mon, z_mon]);
  let two_plus_z = PolyOwning::new(~[2.,1.], ~[one_mon, z_mon]);

  assert_eq!(mul_polys(&one_plus_z, &one_plus_z), 
             PolyOwning::new(~[1.,1.,1.,1.], ~[one_mon, z_mon, z_mon, z_mon*z_mon]));
  assert_eq!(mul_polys(&one_plus_z, &two_plus_z),
             PolyOwning::new(~[2.,1.,2.,1.], ~[one_mon, z_mon, z_mon, z_mon*z_mon]));
}


#[test]
fn test_3d_owned_canonform() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,1.,1.,1.], ~[z_mon*z_mon, one_mon, z_mon, z_mon]).canonical_form(),
             PolyOwning::new(~[1.,2.,1.], ~[one_mon, z_mon, z_mon*z_mon]));

}

#[test]
fn test_3d_owned_tostr() {
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let y_mon = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };

  assert_eq!(PolyOwning::new(~[1.,2.,3.,4.], ~[z_mon*z_mon, z_mon*y_mon, z_mon, z_mon]).to_str(),
             ~"1 x^0y^0z^2 + 2 x^0y^1z^1 + 3 x^0y^0z^1 + 4 x^0y^0z^1");
}



#[test]
#[should_fail]
fn test_4d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon4d] = ~[Monomial::one()];
  PolyBorrowingMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_4d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon4d] = ~[];
  PolyBorrowingMons::new(no_Rs, no_mons);
}


#[test]
fn test_4d_borrowed_equiv() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, t_mon];
  let one_t_one_t = ~[one_mon, t_mon, one_mon, t_mon];

  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyBorrowingMons::new(~[1.,2.,1.,3.], one_t_one_t)));
  assert!(&PolyBorrowingMons::new(~[2.,5.], mons)
   .equiv(&PolyOwning::new(~[1.,2.,1.,3.], one_t_one_t)));
}  


#[test]
fn test_4d_borrowed_scaling() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, t_mon];
  assert_eq!(PolyBorrowingMons::new(~[1.,2.], mons).scaled(2.),
             PolyBorrowingMons::new(~[2.,4.], mons));
}


#[test]
fn test_4d_borrowed_lcomb() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let t_mon_singleton = ~[t_mon];
  let one = PolyBorrowingMons::new(~[1.], one_mon_singleton);
  let t = PolyBorrowingMons::new(~[1.], t_mon_singleton);

  assert_eq!(PolyOwning::from_polys_lcomb([(-2.,&one), (3.,&t), (-2.5,&t), (1.,&one)]), poly([(-1.,one_mon),(0.5,t_mon)]));
}


#[test]
#[should_fail]
fn test_4d_owned_constr_mismatched() {
  let one_mon: Mon4d = Monomial::one();
  PolyOwning::new(~[1.,1.], ~[one_mon]);
}

#[test]
#[should_fail]
fn test_4d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon4d] = ~[];
  PolyOwning::new(no_Rs, no_mons);
}


#[test]
fn test_4d_owned_equiv() {
  let one: Mon4d = Monomial::one();
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert!(poly([(1.,one),(2.,t), (1.,one), (3.,t)])
            .equiv(&poly([(2.,one),(5.,t)])));

  let mons = ~[one, t];
  assert!(poly([(1.,one),(2.,t), (1.,one), (3.,t)])
            .equiv(&PolyBorrowingMons::new(~[2.,5.], mons)));
}  

#[test]
fn test_4d_owned_scaling() {
  let one: Mon4d = Monomial::one();
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,2.], ~[one, t]).scaled(2.),
             PolyOwning::new(~[2.,4.], ~[one, t]))
}

#[test]
fn test_4d_owned_lcomb() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one = PolyOwning::new(~[1.], ~[one_mon]);
  let t = PolyOwning::new(~[1.], ~[t_mon]);

  assert_eq!(PolyOwning::from_polys_lcomb([(-2.,&one), (3.,&t), (-2.5,&t), (1.,&one)]), poly([(-1.,one_mon),(0.5,t_mon)]));
}

#[test]
fn test_4d_owned_multiplication() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one_plus_t = PolyOwning::new(~[1.,1.], ~[one_mon, t_mon]);
  let two_plus_t = PolyOwning::new(~[2.,1.], ~[one_mon, t_mon]);

  assert_eq!(mul_polys(&one_plus_t, &one_plus_t), 
             PolyOwning::new(~[1.,1.,1.,1.], ~[one_mon, t_mon, t_mon, t_mon*t_mon]));
  assert_eq!(mul_polys(&one_plus_t, &two_plus_t),
             PolyOwning::new(~[2.,1.,2.,1.], ~[one_mon, t_mon, t_mon, t_mon*t_mon]));
}


#[test]
fn test_4d_owned_canonform() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert_eq!(PolyOwning::new(~[1.,1.,1.,1.], ~[t_mon*t_mon, one_mon, t_mon, t_mon]).canonical_form(),
             PolyOwning::new(~[1.,2.,1.], ~[one_mon, t_mon, t_mon*t_mon]));
}

#[test]
fn test_4d_owned_tostr() {
  let x_mon = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y_mon = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };

  assert_eq!(PolyOwning::new(~[1.,2.,3.,4.], ~[x_mon*x_mon, x_mon*y_mon, x_mon, x_mon]).to_str(),
             ~"1 x1^2x2^0x3^0x4^0 + 2 x1^1x2^1x3^0x4^0 + 3 x1^1x2^0x3^0x4^0 + 4 x1^1x2^0x3^0x4^0");
}


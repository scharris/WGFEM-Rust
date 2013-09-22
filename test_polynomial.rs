use common::*;
use monomial::{Monomial, Mon1d, Mon2d, Mon3d, Mon4d};
use polynomial::*;

#[test]
#[should_fail]
fn test_1d_owned_constr_mismatched() {
  let one_mon: Mon1d = Monomial::one();
  PolyWithOwnedMons::new(~[1.,1.], ~[one_mon]);
}

#[test]
#[should_fail]
fn test_1d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon1d] = ~[];
  PolyWithOwnedMons::new(no_Rs, no_mons);
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
  PolyWithBorrowedMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_1d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon1d] = ~[];
  PolyWithBorrowedMons::new(no_Rs, no_mons);
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


#[test]
#[should_fail]
fn test_2d_owned_constr_mismatched() {
  let one_mon: Mon2d = Monomial::one();
  PolyWithOwnedMons::new(~[1.,1.], ~[one_mon]);
}


#[test]
#[should_fail]
fn test_2d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon2d] = ~[];
  PolyWithOwnedMons::new(no_Rs, no_mons);
}


#[test]
fn test_2d_owned_equiv() {
  let one: Mon2d = Monomial::one();
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  assert!(poly(~[(1.,one),(2.,y), (1.,one), (3.,y)])
            .equiv(&poly(~[(2.,one),(5.,y)])));

  let mons = ~[one, y];
  assert!(poly(~[(1.,one),(2.,y), (1.,one), (3.,y)])
            .equiv(&PolyWithBorrowedMons::new(~[2.,5.], mons)));
}  

#[test]
fn test_2d_owned_scaling() {
  let one: Mon2d = Monomial::one();
  let y = Mon2d { exps: [Deg(1), Deg(0)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.], ~[one, y]).scaled(2.),
             PolyWithOwnedMons::new(~[2.,4.], ~[one, y]))
  // scaling by mutation
  let mut p = PolyWithOwnedMons::new(~[1.,2.], ~[one, y]);
  p.scale(3.);
  assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,6.], ~[one, y])));
}  

#[test]
fn test_2d_owned_addition() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let y = PolyWithOwnedMons::new(~[1.], ~[y_mon]);

  assert_eq!(one + y, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,y_mon]));
  assert_eq!(one + one + y, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, one_mon, y_mon]));
  assert_eq!(one + y + y, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, y_mon, y_mon]));
}

#[test]
fn test_2d_owned_multiplication() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one_plus_y = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, y_mon]);
  let two_plus_y = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, y_mon]);

  assert_eq!(one_plus_y * one_plus_y, 
             PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[one_mon, y_mon, y_mon, y_mon*y_mon]));
  assert_eq!(one_plus_y * two_plus_y,
             PolyWithOwnedMons::new(~[2.,1.,2.,1.], ~[one_mon, y_mon, y_mon, y_mon*y_mon]));
}


#[test]
fn test_2d_owned_canonform() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[y_mon*y_mon, one_mon, y_mon, y_mon]).canonical_form(),
             PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon, y_mon, y_mon*y_mon]));

}

#[test]
fn test_2d_owned_tostr() {
  let x_mon = Mon2d { exps: [Deg(1), Deg(0)] };
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.,3.,4.], ~[x_mon*x_mon, x_mon*y_mon, x_mon, x_mon]).to_str(),
             ~"1 x^2y^0 + 2 x^1y^1 + 3 x^1y^0 + 4 x^1y^0");
}



#[test]
#[should_fail]
fn test_2d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon2d] = ~[Monomial::one()];
  PolyWithBorrowedMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_2d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon2d] = ~[];
  PolyWithBorrowedMons::new(no_Rs, no_mons);
}


#[test]
fn test_2d_borrowed_equiv() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let mons = ~[one_mon, y_mon];
  let one_y_one_y = ~[one_mon, y_mon, one_mon, y_mon];

  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithBorrowedMons::new(~[1.,2.,1.,3.], one_y_one_y)));
  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.,3.], one_y_one_y)));
}  


#[test]
fn test_2d_borrowed_scaling() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let mons = ~[one_mon, y_mon];
  assert_eq!(PolyWithBorrowedMons::new(~[1.,2.], mons).scaled(2.),
             PolyWithBorrowedMons::new(~[2.,4.], mons));
  // scaling by mutation
  let mut p = PolyWithBorrowedMons::new(~[1.,2.], mons);
  p.scale(3.);
  assert!(p.equiv(&PolyWithBorrowedMons::new(~[3.,6.], mons)));
}  


#[test]
fn test_2d_borrowed_addition() {
  let one_mon: Mon2d = Monomial::one();
  let y_mon = Mon2d { exps: [Deg(0), Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let y_mon_singleton = ~[y_mon];
  let one = PolyWithBorrowedMons::new(~[1.], one_mon_singleton);
  let y = PolyWithBorrowedMons::new(~[1.], y_mon_singleton);

  assert_eq!(one + y, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, y_mon]));
}


#[test]
#[should_fail]
fn test_3d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon3d] = ~[Monomial::one()];
  PolyWithBorrowedMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_3d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon3d] = ~[];
  PolyWithBorrowedMons::new(no_Rs, no_mons);
}

#[test]
fn test_3d_borrowed_equiv() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, z_mon];
  let one_z_one_z = ~[one_mon, z_mon, one_mon, z_mon];

  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithBorrowedMons::new(~[1.,2.,1.,3.], one_z_one_z)));
  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.,3.], one_z_one_z)));
}  


#[test]
fn test_3d_borrowed_scaling() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, z_mon];
  assert_eq!(PolyWithBorrowedMons::new(~[1.,2.], mons).scaled(2.),
             PolyWithBorrowedMons::new(~[2.,4.], mons));
  // scaling by mutation
  let mut p = PolyWithBorrowedMons::new(~[1.,2.], mons);
  p.scale(3.);
  assert!(p.equiv(&PolyWithBorrowedMons::new(~[3.,6.], mons)));
}  


#[test]
fn test_3d_borrowed_addition() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let z_mon_singleton = ~[z_mon];
  let one = PolyWithBorrowedMons::new(~[1.], one_mon_singleton);
  let z = PolyWithBorrowedMons::new(~[1.], z_mon_singleton);

  assert_eq!(one + z, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, z_mon]));
}

#[test]
#[should_fail]
fn test_3d_owned_constr_mismatched() {
  let one_mon: Mon3d = Monomial::one();
  PolyWithOwnedMons::new(~[1.,1.], ~[one_mon]);
}

#[test]
#[should_fail]
fn test_3d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon3d] = ~[];
  PolyWithOwnedMons::new(no_Rs, no_mons);
}


#[test]
fn test_3d_owned_equiv() {
  let one: Mon3d = Monomial::one();
  let z = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };

  assert!(poly(~[(1.,one),(2.,z), (1.,one), (3.,z)])
            .equiv(&poly(~[(2.,one),(5.,z)])));

  let mons = ~[one, z];
  assert!(poly(~[(1.,one),(2.,z), (1.,one), (3.,z)])
            .equiv(&PolyWithBorrowedMons::new(~[2.,5.], mons)));
}  

#[test]
fn test_3d_owned_scaling() {
  let one: Mon3d = Monomial::one();
  let z = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.], ~[one, z]).scaled(2.),
             PolyWithOwnedMons::new(~[2.,4.], ~[one, z]))
  // scaling by mutation
  let mut p = PolyWithOwnedMons::new(~[1.,2.], ~[one, z]);
  p.scale(3.);
  assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,6.], ~[one, z])));
}  

#[test]
fn test_3d_owned_addition() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let z = PolyWithOwnedMons::new(~[1.], ~[z_mon]);

  assert_eq!(one + z, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,z_mon]));
  assert_eq!(one + one + z, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, one_mon, z_mon]));
  assert_eq!(one + z + z, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, z_mon, z_mon]));
}

#[test]
fn test_3d_owned_multiplication() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_plus_z = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, z_mon]);
  let two_plus_z = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, z_mon]);

  assert_eq!(one_plus_z * one_plus_z, 
             PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[one_mon, z_mon, z_mon, z_mon*z_mon]));
  assert_eq!(one_plus_z * two_plus_z,
             PolyWithOwnedMons::new(~[2.,1.,2.,1.], ~[one_mon, z_mon, z_mon, z_mon*z_mon]));
}


#[test]
fn test_3d_owned_canonform() {
  let one_mon: Mon3d = Monomial::one();
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[z_mon*z_mon, one_mon, z_mon, z_mon]).canonical_form(),
             PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon, z_mon, z_mon*z_mon]));

}

#[test]
fn test_3d_owned_tostr() {
  let z_mon = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let y_mon = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.,3.,4.], ~[z_mon*z_mon, z_mon*y_mon, z_mon, z_mon]).to_str(),
             ~"1 x^0y^0z^2 + 2 x^0y^1z^1 + 3 x^0y^0z^1 + 4 x^0y^0z^1");
}



#[test]
#[should_fail]
fn test_4d_borrowed_constr_mismatched() {
  let one_mon_singleton: ~[Mon4d] = ~[Monomial::one()];
  PolyWithBorrowedMons::new(~[1.,1.], one_mon_singleton);
}

#[test]
#[should_fail]
fn test_4d_borrowed_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon4d] = ~[];
  PolyWithBorrowedMons::new(no_Rs, no_mons);
}


#[test]
fn test_4d_borrowed_equiv() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, t_mon];
  let one_t_one_t = ~[one_mon, t_mon, one_mon, t_mon];

  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithBorrowedMons::new(~[1.,2.,1.,3.], one_t_one_t)));
  assert!(&PolyWithBorrowedMons::new(~[2.,5.], mons)
   .equiv(&PolyWithOwnedMons::new(~[1.,2.,1.,3.], one_t_one_t)));
}  


#[test]
fn test_4d_borrowed_scaling() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let mons = ~[one_mon, t_mon];
  assert_eq!(PolyWithBorrowedMons::new(~[1.,2.], mons).scaled(2.),
             PolyWithBorrowedMons::new(~[2.,4.], mons));
  // scaling by mutation
  let mut p = PolyWithBorrowedMons::new(~[1.,2.], mons);
  p.scale(3.);
  assert!(p.equiv(&PolyWithBorrowedMons::new(~[3.,6.], mons)));
}  


#[test]
fn test_4d_borrowed_addition() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one_mon_singleton = ~[one_mon];
  let t_mon_singleton = ~[t_mon];
  let one = PolyWithBorrowedMons::new(~[1.], one_mon_singleton);
  let t = PolyWithBorrowedMons::new(~[1.], t_mon_singleton);

  assert_eq!(one + t, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, t_mon]));
}


#[test]
#[should_fail]
fn test_4d_owned_constr_mismatched() {
  let one_mon: Mon4d = Monomial::one();
  PolyWithOwnedMons::new(~[1.,1.], ~[one_mon]);
}

#[test]
#[should_fail]
fn test_4d_owned_constr_zero_size() {
  let no_Rs: ~[R] = ~[];
  let no_mons: ~[Mon4d] = ~[];
  PolyWithOwnedMons::new(no_Rs, no_mons);
}


#[test]
fn test_4d_owned_equiv() {
  let one: Mon4d = Monomial::one();
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert!(poly(~[(1.,one),(2.,t), (1.,one), (3.,t)])
            .equiv(&poly(~[(2.,one),(5.,t)])));

  let mons = ~[one, t];
  assert!(poly(~[(1.,one),(2.,t), (1.,one), (3.,t)])
            .equiv(&PolyWithBorrowedMons::new(~[2.,5.], mons)));
}  

#[test]
fn test_4d_owned_scaling() {
  let one: Mon4d = Monomial::one();
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.], ~[one, t]).scaled(2.),
             PolyWithOwnedMons::new(~[2.,4.], ~[one, t]))
  // scaling by mutation
  let mut p = PolyWithOwnedMons::new(~[1.,2.], ~[one, t]);
  p.scale(3.);
  assert!(p.equiv(&PolyWithOwnedMons::new(~[3.,6.], ~[one, t])));
}  

#[test]
fn test_4d_owned_addition() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one = PolyWithOwnedMons::new(~[1.], ~[one_mon]);
  let t = PolyWithOwnedMons::new(~[1.], ~[t_mon]);

  assert_eq!(one + t, PolyWithOwnedMons::new(~[1.,1.], ~[one_mon,t_mon]));
  assert_eq!(one + one + t, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, one_mon, t_mon]));
  assert_eq!(one + t + t, PolyWithOwnedMons::new(~[1.,1.,1.], ~[one_mon, t_mon, t_mon]));
}

#[test]
fn test_4d_owned_multiplication() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one_plus_t = PolyWithOwnedMons::new(~[1.,1.], ~[one_mon, t_mon]);
  let two_plus_t = PolyWithOwnedMons::new(~[2.,1.], ~[one_mon, t_mon]);

  assert_eq!(one_plus_t * one_plus_t, 
             PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[one_mon, t_mon, t_mon, t_mon*t_mon]));
  assert_eq!(one_plus_t * two_plus_t,
             PolyWithOwnedMons::new(~[2.,1.,2.,1.], ~[one_mon, t_mon, t_mon, t_mon*t_mon]));
}


#[test]
fn test_4d_owned_canonform() {
  let one_mon: Mon4d = Monomial::one();
  let t_mon = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,1.,1.,1.], ~[t_mon*t_mon, one_mon, t_mon, t_mon]).canonical_form(),
             PolyWithOwnedMons::new(~[1.,2.,1.], ~[one_mon, t_mon, t_mon*t_mon]));
}

#[test]
fn test_4d_owned_tostr() {
  let x_mon = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y_mon = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };

  assert_eq!(PolyWithOwnedMons::new(~[1.,2.,3.,4.], ~[x_mon*x_mon, x_mon*y_mon, x_mon, x_mon]).to_str(),
             ~"1 x1^2x2^0x3^0x4^0 + 2 x1^1x2^1x3^0x4^0 + 3 x1^1x2^0x3^0x4^0 + 4 x1^1x2^0x3^0x4^0");
}


extern mod extra;
use extra::treemap::TreeMap;
use std::vec;

type R = float;
type Deg = u8;
fn deg(n: int) -> Deg { n as Deg }

type Dim = u8;
fn dim(d: int) -> Dim { d as Dim }


#[deriving(Eq, IterBytes, TotalOrd, TotalEq, Clone, Ord)]
pub struct Monomial {
  exps: ~[Deg]
}


impl Monomial {

  fn with_exps(exps: ~[Deg]) -> Monomial {
    if exps.len() == 0 { fail!("Monomial requires non-empty array of exponents."); }
    else {
      Monomial {exps: exps}
    }
  }

  fn domain_dim(&self) -> Dim {
    self.exps.len() as Dim
  }

  fn one(d: Dim) -> Monomial {
    let exps = vec::from_elem(d as uint, 0 as Deg);
    Monomial::with_exps(exps)
  }
}

impl ToStr for Monomial {
  fn to_str(&self) -> ~str {
    match self.domain_dim() {
      1 => { fmt!("x^%?", self.exps[0]) }
      2 => { fmt!("x^%?y^%?", self.exps[0], self.exps[1]) }
      3 => { fmt!("x^%?y^%?z^%?", self.exps[0], self.exps[1], self.exps[2]) }
      _ => { /* join(map(i->"x$i^$(int(m.exps[i]))", 1:ddim), " ") */ ~"TODO" }
    }
  }
}

#[deriving(Eq, IterBytes, Clone, Ord)]
pub struct Polynomial {
  mons: ~[Monomial],
  coefs: ~[R]
}



impl Polynomial {

  fn with_mons_and_coefs(mons: ~[Monomial], coefs: ~[R]) -> Polynomial {
    if mons.len() != coefs.len() { fail!("size of coefficent and monomial arrays must match"); }
    else if mons.len() == 0 { fail!("polynomial requires one or more terms"); }
    else {
      Polynomial { mons: mons, coefs: coefs }
    }
  }

  fn canonical_form(&self) -> Polynomial {

    let combined_coefs_by_mon = {
      let mut mon_to_coef:TreeMap<&Monomial,R> = TreeMap::new();
      for i in range(0, self.mons.len()) {
        let (mon, mon_coef) = (&self.mons[i], self.coefs[i]);
        if mon_coef != 0 as R {
          let updated = match mon_to_coef.find_mut(&mon) {
            Some(coef) => { *coef += mon_coef; true }, None => false
          };
          if !updated {
            mon_to_coef.insert(mon, mon_coef);
          }
        }
      }
      mon_to_coef
    };

    let mut mons:~[Monomial] = vec::with_capacity(combined_coefs_by_mon.len());
    let mut coefs:~[R] = vec::with_capacity(combined_coefs_by_mon.len());
    for (&mon, &coef) in combined_coefs_by_mon.iter() {
      if coef != 0 as R {
        mons.push(mon.clone());
        coefs.push(coef);
      }
    }

    if mons.len() != 0 {
      Polynomial::with_mons_and_coefs(mons, coefs)
    }
    else {
      Polynomial::zero(self.domain_dim())
    }
  }

  fn domain_dim(&self) -> Dim {
    self.mons[0].domain_dim()
  }

  fn zero(d: Dim) -> Polynomial {
    Polynomial::with_mons_and_coefs(~[Monomial::one(d)], ~[0 as R])
  }
}

impl ToStr for Polynomial {
  fn to_str(&self) -> ~str {
    // join(map(i->"$(p.coefs[i]) $(p.mons[i])", 1:length(p.mons)), " + ")
    fmt!("Polynomial %?", *self)
  }
}

fn main() {
  let x = Monomial::with_exps(~[deg(1), deg(0)]);
  let y = Monomial::with_exps(~[deg(0), deg(1)]);
  let xy = Monomial::with_exps(~[deg(1), deg(1)]);
  let mons = ~[x.clone(), y.clone(), xy.clone(), x.clone(), xy.clone()];
  let p = Polynomial::with_mons_and_coefs(mons, ~[1.0, -1.0, 2.0, -1.0, 2.0]);
  let cfp = p.canonical_form();
  println(fmt!("Canonical form of %s is %s", p.to_str(), cfp.to_str()));
}

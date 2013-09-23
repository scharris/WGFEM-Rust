use common::*;
use std::vec;

// row major dense matrix type
struct DenseMatrix {
  num_rows: uint,
  num_cols: uint,
  data: ~[R]
}


impl DenseMatrix {
  
  pub fn symmetric_from_lower_triangular_fn(num_rows: uint, num_cols: uint, f: &fn(row:uint, col:uint) -> R) -> DenseMatrix {
    let mut data = vec::from_elem(num_rows * num_cols, 0 as R);
    for r in range(0, num_rows) {
      for c in range(0, r) {
        let f_val = f(r,c);
        data[r * num_cols + c] = f_val;
        data[c * num_cols + r] = f_val;
      }
      data[r * num_cols + r] = f(r,r); 
    }
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: data
    }
  }

  pub fn from_fn(num_rows: uint, num_cols: uint, f: &fn(row:uint, col:uint) -> R) -> DenseMatrix {
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: vec::from_fn(num_rows * num_cols, |i| {
        let r = i / num_cols;
        let c = i % num_cols;
        f(r,c)
      })
    }
  }

  /*
  pub fn from_elem(num_rows: uint, num_cols: uint, elem: R) -> DenseMatrix {
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: vec::from_elem(num_rows * num_cols, elem)
    }
  }
  */


  #[inline]
  pub fn get(&self, r: uint, c: uint) -> R {
    self.data[r * self.num_cols + c]
  }
  
  #[inline]
  pub fn set(&mut self, r: uint, c: uint, value: R) {
    self.data[r * self.num_cols + c] = value;
  }
}


#[test]
fn test_constr_from_fn() {
  let m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.data, ~[0.,1.,2.,1000.,1001.,1002.,2000.,2001.,2002.]);
}

#[test]
fn test_constr_symmetric_from_lower_triangular_fn() {
  let m = DenseMatrix::symmetric_from_lower_triangular_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.data, ~[0.,1000.,2000.,1000.,1001.,2001.,2000.,2001.,2002.]);
}

#[test]
fn test_get() {
  let m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,0), 1000.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,0), 2000.);
  assert_eq!(m.get(2,1), 2001.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
fn test_set() {
  let mut m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.get(2,1), 2001.);
  m.set(2,1, 2001.5);
  assert_eq!(m.get(2,1), 2001.5);
}


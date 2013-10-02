use common::*;
use std::vec;

/// Column major dense matrix type.
pub struct DenseMatrix {
  num_rows: uint,
  num_cols: uint,
  data: ~[R]
}


impl DenseMatrix {
 
  pub fn from_fn(num_rows: uint, num_cols: uint, f: &fn(row:uint, col:uint) -> R) -> DenseMatrix {
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: vec::from_fn(num_rows * num_cols, |i| {
        let c = i / num_rows;
        let r = i % num_rows;
        f(r,c)
      })
    }
  }
  
  pub fn from_elem(num_rows: uint, num_cols: uint, elem: R) -> DenseMatrix {
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: vec::from_elem(num_rows * num_cols, elem)
    }
  }

  pub fn upper_triangular_from_fn(num_rows: uint, num_cols: uint, f: &fn(row:uint, col:uint) -> R) -> DenseMatrix {
    let mut data = vec::from_elem(num_rows * num_cols, 0 as R);
    for c in range(0, num_cols) {
      for r in range(0, c) {
        let f_val = f(r,c);
        data[c * num_rows + r] = f_val;
      }
      data[c * num_cols + c] = f(c,c); 
    }
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: data
    }
  }

  #[inline]
  pub fn get(&self, r: uint, c: uint) -> R {
    self.data[c * self.num_rows + r]
  }
  
  #[inline]
  pub fn set(&mut self, r: uint, c: uint, value: R) {
    self.data[c * self.num_rows + r] = value;
  }

  pub fn col_maj_data_copy(&self) -> ~[R] {
    self.data.clone()
  }

/*
  #[inline]
  pub fn copy_data_into(&self, buf: &mut [R]) -> () {
    use std::vec::MutableCloneableVector;
    let ncopied = buf.copy_from(self.data);
    assert!(ncopied == self.data.len());
  }
*/
}


#[test]
fn test_constr_from_fn() {
  let m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.data, ~[0.,1000.,2000.,1.,1001.,2001.,2.,1002.,2002.]);
}

fn test_constr_from_elem() {
  let m = DenseMatrix::from_elem(3,3, 10.);
  assert_eq!(m.data, ~[10.,10.,10.,10.,10.,10.,10.,10.,10.]);
}


#[test]
fn test_constr_upper_triangular_from_fn() {
  let m = DenseMatrix::upper_triangular_from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.data, ~[0.,0.,0.,1.,1001.,0.,2.,1002.,2002.]);
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


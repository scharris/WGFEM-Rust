use common::*;
use std::vec;

// row major dense matrix type
struct DenseMatrix {
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
        let r = i / num_cols;
        let c = i % num_cols;
        f(r,c)
      })
    }
  }

  /*pub fn from_elem(num_rows: uint, num_cols: uint, elem: R) -> DenseMatrix {
    DenseMatrix {
      num_rows: num_rows,
      num_cols: num_cols,
      data: vec::from_elem(num_rows * num_cols, elem)
    }
  }*/


  #[inline]
  pub fn get(&self, r: uint, c: uint) -> R {
    self.data[r * self.num_cols + c]
  }
  
  #[inline]
  pub fn set(&mut self, r: uint, c: uint, value: R) {
    self.data[r * self.num_cols + c] = value;
  }
}


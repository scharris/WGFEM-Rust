use common::{R};
use lapack::{lapack_int};

use std::vec;
use std::vec::raw;
use std::cast;
use std::ptr;

/// Sparse matrix type, with compressed sparse row storage, 3-array variation (CSR3).
/// Values must be pushed into the matrix in increasing order of their (row, column)
/// pairs, with the row being most significant, and with each row being represented by
/// at least one pushed value (which may be 0).
pub struct SparseMatrix {
  values: ~[R],
  value_cols: ~[lapack_int],
  row_begins: ~[lapack_int],
}

impl SparseMatrix {

  pub fn new_with_capacities(values_capacity: uint, rows_capacity: uint) -> SparseMatrix {
    SparseMatrix {
      values: vec::with_capacity(values_capacity),
      value_cols: vec::with_capacity(values_capacity),
      row_begins: vec::with_capacity(rows_capacity+1),
    }
  }

  #[inline]
  pub fn push(&mut self, r: uint, c: uint, val: R) {
    match r {
      last_row if last_row == self.row_begins.len()-1 => {
        if c <= self.value_cols[self.value_cols.len()-1] as uint {
          fail!("Column numbers must be pushed in strictly increasing order within sparse matrix rows.")
        }
      }
      new_row if new_row == self.row_begins.len() => {
        self.row_begins.push(self.values.len() as lapack_int);
      }
      _ => fail!("Push to sparse matrix must be to last existing row or next non-existing row")
    }
    self.values.push(val);
    self.value_cols.push(c as lapack_int);
  }

  pub fn num_rows(&self) -> uint {
    self.row_begins.len()
  }
  
  pub fn num_values(&self) -> uint {
    self.values.len()
  }
  
  pub fn get(&self, r: uint, c: uint) -> R {
    let next_row_begin = if r == self.row_begins.len()-1 { self.values.len()} else { self.row_begins[r+1] as uint };
    for i in range(self.row_begins[r] as uint, next_row_begin) {
      match (c as lapack_int).cmp(&self.value_cols[i]) {
        Less =>  { return 0 as R; }
        Equal => { return self.values[i]; }
        Greater => {} // keep looking 
      }
    }
    0 as R
  }

  /// Returns the 3-array variant (a, ia, ja) of the Compressed Sparse Row format as pointers to the values array a,
  /// the row beginning indexes ia into the values array, and the corresponding column numbers ja of the values.
  /// The row_begins vector must have a capacity of at least one greater than its length, otherwise an error
  /// is generated. This allows an extra "cap" entry to be written past the proper row beginning index values
  /// as required by lapack, in the reserved capacity part of the row_begins buffer.
  pub unsafe fn csr3_ptrs(&self) -> (*R, *lapack_int, *lapack_int) {
    if self.row_begins.capacity() <= self.row_begins.len() {
      fail!("Need one capacity element beyond the length of the row_begins vector for conversion to raw pointers.");
    }
    else {
      // Cap the row_begins buffer with the number of values as trailing "ghost" row start index, as required by lapack. 
      let row_begins_ptr = raw::to_mut_ptr(cast::transmute_mut(self).row_begins);
      *ptr::mut_offset(row_begins_ptr, self.row_begins.len() as int) = self.values.len() as lapack_int;
      
      (raw::to_ptr(self.values),
       raw::to_ptr(self.row_begins),
       raw::to_ptr(self.value_cols))
    }
  }

}


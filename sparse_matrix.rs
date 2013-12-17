use common::{R};
use la;
use la::lapack_int;

use extra::c_vec::CVec;
use std::cast::transmute;
use std::ptr;
use std::libc::{c_ulong};

/// Sparse matrix type, with compressed sparse row storage, 3-array variation (CSR3).
/// Values must be pushed into the matrix in increasing order of their (row, column)
/// pairs, with the row being most significant, and with each row being represented by
/// at least one pushed value (which may be 0).
pub struct SparseMatrix {

  priv values: CVec<R>,
  priv value_cols: CVec<lapack_int>,
  priv row_first_value_ixs: CVec<lapack_int>,

  priv num_values: uint,
  priv num_rows: uint,

  priv matrix_type: MatrixType,
}

pub enum MatrixType {
  Symmetric,             // symmetric with values in upper triangle
  StructurallySymmetric, // structurally symmetric, m_{i,j} present iff m_{j,i} present.
  General,
}

impl SparseMatrix {

  pub fn new_with_capacities(values_capacity: uint, rows_capacity: uint, mtype: MatrixType) -> SparseMatrix {
    let (values, value_cols, row_first_value_ixs) = unsafe {
      (CVec::new(la::alloc_doubles(values_capacity as c_ulong), values_capacity),
       CVec::new(la::alloc_ints(values_capacity as c_ulong), values_capacity),
       CVec::new(la::alloc_ints((rows_capacity+1u) as c_ulong), rows_capacity)) // alloc extra element for cap value
    };
    SparseMatrix {
      values: values,
      value_cols: value_cols,
      row_first_value_ixs: row_first_value_ixs,
      num_values: 0u,
      num_rows: 0u,
      matrix_type: mtype
    }
  }

  #[inline]
  pub fn push(&mut self, r: uint, c: uint, val: R) {
    match r {
      // If continuing on the same row, the column number should be greater than the last.
      last_row if last_row == self.num_rows-1 => {
        if c <= *self.value_cols.get(self.num_values-1) as uint {
          fail!("Columns must be added in strictly increasing order within sparse matrix rows.")
        }
      }
      // Else the new row must be the next row sequentially. Add the initial value index for the row.
      new_row if new_row == self.num_rows => {
        *self.row_first_value_ixs.get_mut(new_row) = self.num_values as lapack_int;
        self.num_rows += 1;
      }
      _ => fail!("Push to sparse matrix must be to last existing row or next non-existing row")
    }
   
    // Add the value and its column index.
    *self.values.get_mut(self.num_values) = val;
    *self.value_cols.get_mut(self.num_values) = c as lapack_int;
    self.num_values += 1;
  }

  pub fn num_rows(&self) -> uint {
    self.num_rows
  }
  
  pub fn num_values(&self) -> uint {
    self.num_values
  }
  
  pub fn matrix_type(&self) -> MatrixType { self.matrix_type }
  
  pub fn get(&self, r: uint, c: uint) -> R {
    if r >= self.num_rows { fail!("Row index out of range.") }
    let first_val_ix = *self.row_first_value_ixs.get(r) as uint;
    let next_row_begin = if r == self.num_rows-1 { self.num_values } else { *self.row_first_value_ixs.get(r+1) as uint };
    for i in range(first_val_ix, next_row_begin) {
      match (c as lapack_int).cmp(self.value_cols.get(i)) {
        Less =>  { return 0 as R; }
        Equal => { return *self.values.get(i); }
        Greater => {} // keep looking 
      }
    }
    0 as R
  }

  pub fn debug_print(&self) {
    unsafe {
      for r in range(0, self.num_rows) {
        let first_val_ix = *self.row_first_value_ixs.get(r) as uint;
        let next_row_begin = if r == self.num_rows-1 { self.num_values } else { *self.row_first_value_ixs.get(r+1) as uint };
        for i in range(first_val_ix, next_row_begin) {
          let c = *self.value_cols.get(i);
          println!("{}\t{}\t{:.10f}", r+1, c+1, *self.values.get(i));
        }
      }
    }
  }

  /// Returns the 3-array variant (a, ia, ja) of the Compressed Sparse Row format as pointers to the values array a,
  /// the row beginning indexes ia into the values array, and the corresponding column numbers ja of the values.
  /// The row_first_value_ixs vector must have a capacity of at least one greater than its length, otherwise an error
  /// is generated. This allows an extra "cap" entry to be written past the proper row beginning index values
  /// as required by lapack, in the reserved capacity part of the row_first_value_ixs buffer.
  pub unsafe fn csr3_ptrs(&self) -> (*R, *lapack_int, *lapack_int) {
    // Cap the row_first_value_ixs buffer with the number of values as required by some solvers (other solvers unaffected).
    // Extra storage for this item was allocated and is gauranteed to still be unused because it was not made available
    // through the bounds-checked cvec accessors used elsewhere in this implementation.
    *ptr::mut_offset(transmute(self.row_first_value_ixs.get(0)), self.num_rows as int) = self.num_values as lapack_int;

    (transmute(self.values.get(0)),
     transmute(self.row_first_value_ixs.get(0)),
     transmute(self.value_cols.get(0)))
  }

}

#[unsafe_destructor]
impl Drop for SparseMatrix {
  #[inline(never)]
  fn drop(&mut self) {
    unsafe {
      la::free_doubles(transmute(self.values.get(0)));
      la::free_ints(transmute(self.row_first_value_ixs.get(0)));
      la::free_ints(transmute(self.value_cols.get(0)));
    }
  }
}


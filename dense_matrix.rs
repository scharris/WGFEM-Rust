use common::*;
use la;

use std::libc::{c_ulong};
use std::ptr;
use std::iter::{range_inclusive};
use std::cast::transmute;
use extra::c_vec::CVec;

/// Column major dense matrix type.
pub struct DenseMatrix {
  priv data: CVec<R>,
  priv num_rows: uint,
  priv num_cols: uint,
  priv capacity_cols: uint,
}

impl DenseMatrix {

  pub fn from_fn(num_rows: uint, num_cols: uint, f: |row:uint, col:uint| -> R) -> DenseMatrix {
    let n = num_rows * num_cols;
    let mut data = unsafe { alloc_data(n) };
    for i in range(0u, n) {
      let c = i / num_rows;
      let r = i % num_rows;
      unsafe { unsafe_set(&mut data, i, f(r,c)); }
    }
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: num_cols,
    }
  }

  // Create dense matrix of indicated size, with unitialized entries.  This constructor should only be used if the
  // part(s) of the matrix to be used are subsequently initialized, such as via the set() function.
  pub fn of_size(num_rows: uint, num_cols: uint) -> DenseMatrix {
    let n = num_rows * num_cols;
    let data = unsafe { alloc_data(n) };
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: num_cols,
    }
  }
  
  pub fn from_elem(num_rows: uint, num_cols: uint, elem: R) -> DenseMatrix {
    let n = num_rows * num_cols;
    let mut data = unsafe { alloc_data(n) };
    for i in range(0u, n) {
      unsafe { unsafe_set(&mut data, i, elem); }
    }
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: num_cols,
    }
  }
 
  pub fn from_elem_with_cols_capacity(num_rows: uint, num_cols: uint, elem: R, capacity_cols: uint) -> DenseMatrix {
    if capacity_cols < num_cols { fail!("Capacity columns must be greater or equal to number of columns."); }
    let cap = num_rows * capacity_cols;
    let mut data = unsafe { alloc_data(cap) };
    for i in range(0u, cap) {
      unsafe { unsafe_set(&mut data, i, elem); }
    }
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: capacity_cols,
    }
  }
  
  pub fn of_size_with_cols_capacity(num_rows: uint, num_cols: uint, capacity_cols: uint) -> DenseMatrix {
    if capacity_cols < num_cols { fail!("Capacity columns must be greater or equal to number of columns."); }
    let data = unsafe { alloc_data(num_rows * capacity_cols) };
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: capacity_cols,
    }
  }

  // NOTE: This constructor does not initialize the lower triangular values.
  pub fn upper_triangle_from_fn(side_len: uint, f: |row:uint, col:uint| -> R) -> DenseMatrix {
    let n = side_len * side_len;
    let mut data = unsafe { alloc_data(n) };
    for c in range(0, side_len) {
      for r in range_inclusive(0, c) {
        unsafe { unsafe_set(&mut data, c * side_len + r, f(r,c)); }
      }
    }
    DenseMatrix {
      data: data,
      num_rows: side_len,
      num_cols: side_len,
      capacity_cols: side_len,
    }
  }
  
  // NOTE: This constructor does not initialize the upper triangular values.
  pub fn lower_triangle_from_fn(side_len: uint, f: |row:uint, col:uint| -> R) -> DenseMatrix {
    let n = side_len * side_len;
    let mut data = unsafe { alloc_data(n) };
    for r in range(0, side_len) {
      for c in range_inclusive(0, r) {
        unsafe { unsafe_set(&mut data, c * side_len + r, f(r,c)); }
      }
    }
    DenseMatrix {
      data: data,
      num_rows: side_len,
      num_cols: side_len,
      capacity_cols: side_len,
    }
  }

  pub fn from_rows(num_rows: uint, num_cols: uint, elems: &[~[R]]) -> DenseMatrix {
    DenseMatrix::from_fn(num_rows, num_cols, |r,c| elems[r][c])
  }
  
  #[inline(always)]
  pub fn num_rows(&self) -> uint {
    self.num_rows
  }
  
  #[inline(always)]
  pub fn num_cols(&self) -> uint {
    self.num_cols
  }
  
  #[inline(always)]
  pub fn capacity_cols(&self) -> uint {
    self.capacity_cols
  }

  #[inline]
  pub fn get(&self, r: uint, c: uint) -> R {
    if c >= self.num_cols || r >= self.num_rows { fail!("Row or column out of range."); }
    unsafe { unsafe_get(&self.data, c * self.num_rows + r) }
  }
  
  #[inline]
  pub fn set(&mut self, r: uint, c: uint, value: R) {
    if c >= self.num_cols || r >= self.num_rows { fail!("Row or column out of range."); }
    unsafe { unsafe_set(&mut self.data, c * self.num_rows + r, value); }
  }

  #[inline(never)]
  pub fn copy_into(&self, m: &mut DenseMatrix) {
    if self.num_rows != m.num_rows || self.num_cols > m.num_cols {
      fail!("Matrix layouts not compatible for dense matrix copy-into operation.");
    }
    unsafe {
      la::copy_matrix(transmute(self.data.get(0)), self.num_rows as c_ulong, self.num_cols as c_ulong, transmute(m.data.get(0)));
    }
  }
  
  #[inline(never)]
  pub fn copy_upper_triangle_into(&self, m: &mut DenseMatrix) {
    if self.num_rows != m.num_rows || self.num_cols != m.num_cols {
      fail!("Matrix layouts not compatible for dense matrix copy-into operation.");
    }
    unsafe {
      la::copy_upper_triangle(transmute(self.data.get(0)), self.num_rows as c_ulong, self.num_cols as c_ulong, transmute(m.data.get(0)));
    }
  }
  
  #[inline(never)]
  pub fn fill_upper_triangle_from(&mut self, m: &DenseMatrix) {
    if self.num_rows != m.num_rows || m.num_cols > self.num_cols {
      fail!("Matrix layouts not compatible for dense matrix copy-into operation.");
    }
    unsafe {
      la::copy_upper_triangle(transmute(m.data.get(0)), m.num_rows as c_ulong, m.num_cols as c_ulong, transmute(self.data.get(0)));
    }
  }

  #[inline(never)]
  pub fn fill_from_fn(&mut self, f: |row:uint, col:uint| -> R) {
    for r in range(0, self.num_rows) {
      for c in range(0, self.num_cols) {
        unsafe { unsafe_set(&mut self.data, c * self.num_rows + r, f(r,c)); }
      }
    }
  }

  #[inline(always)]
  pub unsafe fn col_maj_data_ptr(&self) -> *R {
    transmute(self.data.get(0))
  }

  #[inline(always)]
  pub unsafe fn mut_col_maj_data_ptr(&self) -> *mut R {
    transmute(self.data.get(0))
  }
  
  pub fn set_num_cols(&mut self, num_cols: uint) {
    if num_cols > self.capacity_cols  {
      fail!("Cannot increase number of columns beyond capacity.");
    }
    self.num_cols = num_cols;
  }

  pub fn print(&self) {
    for i in range(0, self.num_rows) {
      for j in range(0, self.num_cols) {
        print!(" {:6.3f} ", self.get(i,j));
      }
      println!("");
    }
  }
  
} // DenseMatrix impl


#[inline(always)]
unsafe fn unsafe_set(v: &mut CVec<R>, i: uint, value: R) {
  *ptr::mut_offset(transmute(v.get(0)), i as int) = value;
}

#[inline(always)]
unsafe fn unsafe_get(v: &CVec<R>, i: uint) -> R {
  *ptr::mut_offset(transmute(v.get(0)), i as int)
}


#[inline(never)]
pub unsafe fn alloc_data(num_doubles: uint) -> CVec<R> {
  let doubles = la::alloc_doubles(num_doubles as c_ulong);
  CVec::new(doubles, num_doubles)
}

#[unsafe_destructor]
impl Drop for DenseMatrix {
  #[inline(never)]
  fn drop(&mut self) {
    unsafe {
      la::free_doubles(transmute(self.data.get(0)))
    }
  }
}


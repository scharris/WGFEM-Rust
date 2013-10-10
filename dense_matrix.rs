use common::*;
use lapack;

use extra::c_vec;
use extra::c_vec::CVec;
use std::routine::Runnable;
use std::libc::{c_ulong};
use std::ptr;

/// Column major dense matrix type.
pub struct DenseMatrix {
  data: CVec<R>,
  num_rows: uint,
  num_cols: uint,
  capacity_cols: uint,
}

impl DenseMatrix {

  pub fn from_fn(num_rows: uint, num_cols: uint, f: &fn(row:uint, col:uint) -> R) -> DenseMatrix {
    let n = num_rows * num_cols;
    let data = unsafe { alloc_data(n) };
    for i in range(0u, n) {
      let c = i / num_rows;
      let r = i % num_rows;
      c_vec::set(data, i, f(r,c));
    }
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: num_cols,
    }
  }
  
  pub fn from_elem(num_rows: uint, num_cols: uint, elem: R) -> DenseMatrix {
    let n = num_rows * num_cols;
    let data = unsafe { alloc_data(n) };
    for i in range(0u, n) {
      c_vec::set(data, i, elem);
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
    let data = unsafe { alloc_data(cap) };
    for i in range(0u, cap) {
      c_vec::set(data, i, elem);
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

  // NOTE: This constructor does not initialize the lower triangular values: lower triangular values are undefined.
  pub fn with_upper_triangle_from_fn(num_rows: uint, num_cols: uint, f: &fn(row:uint, col:uint) -> R) -> DenseMatrix {
    let n = num_rows * num_cols;
    let data = unsafe { alloc_data(n) };
    for c in range(0, num_cols) {
      for r in range(0, c) {
        c_vec::set(data, c * num_rows + r, f(r,c));
      }
      c_vec::set(data, c * num_cols + c, f(c,c));
    }
    DenseMatrix {
      data: data,
      num_rows: num_rows,
      num_cols: num_cols,
      capacity_cols: num_cols,
    }
  }

  #[inline]
  pub fn get(&self, r: uint, c: uint) -> R {
    if c >= self.num_cols || r >= self.num_rows { fail!("Row or column out of range."); }
    unsafe { *ptr::mut_offset(c_vec::ptr(self.data), (c * self.num_rows + r) as int) }
  }
  
  #[inline]
  pub fn set(&mut self, r: uint, c: uint, value: R) {
    if c >= self.num_cols || r >= self.num_rows { fail!("Row or column out of range."); }
    unsafe { *ptr::mut_offset(c_vec::ptr(self.data), (c * self.num_rows + r) as int) = value; }
  }

  #[fixed_stack_segment]
  #[inline(never)]
  pub fn copy_into(&self, m: &mut DenseMatrix) {
    if self.num_rows != m.num_rows || self.num_cols > m.num_cols {
      fail!("Matrix layouts not compatible for dense matrix copy-into operation.");
    }
    unsafe {
      lapack::copy_matrix(c_vec::ptr(self.data) as *R, self.num_rows as c_ulong, self.num_cols as c_ulong, c_vec::ptr(m.data));
    }
  }
  
  #[fixed_stack_segment]
  #[inline(never)]
  pub fn copy_upper_triangle_into(&self, m: &mut DenseMatrix) {
    if self.num_rows != m.num_rows || self.num_cols > m.num_cols {
      fail!("Matrix layouts not compatible for dense matrix copy-into operation.");
    }
    unsafe {
      lapack::copy_upper_triangle(c_vec::ptr(self.data) as *R, self.num_rows as c_ulong, self.num_cols as c_ulong, c_vec::ptr(m.data));
    }
  }

  #[inline]
  pub unsafe fn col_maj_data_ptr(&self) -> *R {
    c_vec::ptr(self.data) as *R
  }

  #[inline]
  pub unsafe fn mut_col_maj_data_ptr(&self) -> *mut R {
    c_vec::ptr(self.data)
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



#[fixed_stack_segment]
#[inline(never)]
pub unsafe fn alloc_data(num_doubles: uint) -> CVec<R> {
  let doubles = lapack::alloc_doubles(num_doubles as c_ulong);
  let dealloc = ~Deallocator { mem: doubles } as ~Runnable;
  c_vec::c_vec_with_dtor(doubles, num_doubles, dealloc)
}

struct Deallocator {
  mem: *mut R
}

impl Runnable for Deallocator {
  #[fixed_stack_segment]
  fn run(~self) {
    unsafe {
      lapack::free_doubles(self.mem);
    }
  }
}



use std::libc::{c_double, c_ulong, c_int, c_void, malloc, calloc, realloc, free};
use std::cast;

pub type lapack_int = c_int; // Adjust according to whether LP64 or ILP64 libraries are being linked

#[fixed_stack_segment]
#[inline(never)]
pub fn init() {
  unsafe {
    init_allocator(cast::transmute(malloc), cast::transmute(calloc), cast::transmute(realloc), cast::transmute(free));
  }
}

#[link_args = "-Llib/mkl lib/lapack.o -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_core -lmkl_intel_thread -lmkl_core -liomp5"]
extern {

  pub fn init_allocator(malloc_fn: *c_void, calloc_fn: *c_void, realloc_fn: *c_void, free_fn: *c_void);

  pub fn alloc_doubles(num_doubles: c_ulong) -> *mut c_double;

  pub fn free_doubles(mem: *mut c_double);


  pub fn copy_matrix(from_data: *c_double, num_rows: c_ulong, num_cols: c_ulong, to_data: *mut c_double);
  
  pub fn copy_upper_triangle(from_data: *c_double, num_rows: c_ulong, num_cols: c_ulong, to_data: *mut c_double);


  pub fn solve_symmetric_as_col_maj_with_ut_sys(a: *mut c_double,
                                                n: lapack_int,
                                                b: *mut c_double,
                                                nrhs: lapack_int,
                                                ipiv: *mut lapack_int) -> lapack_int;
}


use std::vec;
use std::libc::size_t;
use std::ptr;
use std::libc::{c_uint, c_int, c_void};
use std::cast;

use common::*;

#[fixed_stack_segment]  
#[inline(never)]
pub fn quadrature(f: & &fn(&[R]) -> R, // OK for this to be a closure address
                  min_corner: &[R],
                  max_corner: &[R],
                  rel_err: R, abs_err: R) -> R {
  let (val, status) = unsafe {
    let f_dom_space_dims = min_corner.len() as c_uint;
    let f_range_space_dims = 1 as c_uint;
    let f_pv: *c_void = cast::transmute(f); 
    let integrand_caller_pv: *c_void = cast::transmute(integrand_caller);
    let min_bounds = vec::raw::to_ptr(min_corner);
    let max_bounds = vec::raw::to_ptr(max_corner);
    let max_evals = 0 as size_t;
    let norm_unused = 0u32;
    let mut val = 0 as R;
    let mut err = 0 as R;
    let status = 
      hquadrature(f_range_space_dims,
                  integrand_caller_pv,
                  f_pv,
                  f_dom_space_dims, min_bounds, max_bounds,
                  max_evals, rel_err, abs_err,
                  norm_unused, &mut val, &mut err);
    (val, status)
  };

  if (status != 0) { fail!("quadrature call returned non-zero status"); }
  
  val
}

// This is the integrand callback function called directly by the C integration routine, whose job
// is to calculate the integrand value using f_ptr and set that value in *fval. 
#[inline(never)]
extern fn integrand_caller(ndim: c_uint, x: *R, 
                           f_ptr: *mut (&fn(x:&[R]) -> R), 
                           _: c_uint, fval: *mut R) -> c_int {
  unsafe {
    let f = ptr::read_ptr(f_ptr);
    *fval = vec::raw::buf_as_slice(x, ndim as uint, f);
  }
  0 as c_int
}


// The external C integration routine.
#[link_args = "lib/quadrature.o"]
extern {

  fn hquadrature(fdim: c_uint,
                 f: *c_void, fdata: *c_void,
                 dim: c_uint, xmin: *R, xmax: *R, 
                 maxEval: size_t, reqAbsError: R, reqRelError: R, 
                 norm: u32,
                 val: *mut R, err: *mut R) -> c_int;
}


#[fixed_stack_segment]  
#[test]
fn test_quadrature() {
  let f1 = |_: &[f64]| 2.0;
  let f2 = |x: &[f64]| 2.0*x[0]*x[1];
  let min_bounds = ~[0.,0.];
  let max_bounds = ~[1.,1.];
  assert_eq!(quadrature(&f1, min_bounds, max_bounds, 1e-5, 1e-5), 2.0)
  assert_eq!(quadrature(&f2, min_bounds, max_bounds, 1e-5, 1e-5), 0.5)
  assert_eq!(quadrature(&f1, min_bounds, max_bounds, 1e-5, 1e-5), 2.0)
}


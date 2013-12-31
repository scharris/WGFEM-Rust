use std::vec;
use std::libc::size_t;
use std::ptr;
use std::libc::{c_uint, c_int, c_void};
use std::cast;

use common::*;

#[inline(never)]
pub fn space_adaptive_quadrature(f: & |&[R]| -> R, min_corner: &[R], max_corner: &[R], rel_err: R, abs_err: R) -> R {
  let (val, status) = unsafe {
    let f_dom_space_dims = min_corner.len() as c_uint;
    let f_range_space_dims = 1 as c_uint;
    let f_pv: *c_void = cast::transmute(f); 
    let integrand_caller_pv: *c_void = cast::transmute(h_integrand_caller);
    let min_bounds = min_corner.as_ptr();
    let max_bounds = max_corner.as_ptr();
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

  if (status != 0) { fail!("hquadrature call returned non-zero status"); }
  
  val
}

// Perform Gaussian-Legendre quadrature of f on [a,b]x[c,d] using n weights per axis, for a total of 2n weights and n^2
// evaluation points.
#[inline(never)]
pub fn gaussian_quadrature_2D_rect(n: uint, f: & |x: R, y: R| -> R, a: R, b: R, c: R, d: R) -> R {
  let gq_order = {
    if n <= 10 { 2*n as uint}
    else if n <= 16 { 32 }
    else if n <= 32 { 64 }
    else if n <= 48 { 96 }
    else if n <= 64 { 128 }
    else if n <= 128 { 256 }
    else if n <= 256 { 512 }
    else { fail!("Degree limit exceeded for Gaussian quadrature.") }
  };
  unsafe {
    let f_pv: *c_void = cast::transmute(f); 
    let gq_2D_integrand_caller_pv: *c_void = cast::transmute(gq_2D_integrand_caller);
    gauss_legendre_2D_rect(gq_order as c_int, gq_2D_integrand_caller_pv, f_pv, a, b, c, d)
  }
}

// Perform Gaussian-Legendre quadrature of f on [a,b] using n points.
#[inline(never)]
pub fn gaussian_quadrature(n: uint, f: & |R| -> R, a: R, b: R) -> R {
  let gq_order = {
    if n <= 10 { 2*n as uint}
    else if n <= 16 { 32 }
    else if n <= 32 { 64 }
    else if n <= 48 { 96 }
    else if n <= 64 { 128 }
    else if n <= 128 { 256 }
    else if n <= 256 { 512 }
    else { fail!("Degree limit exceeded for Gaussian quadrature.") }
  };
  unsafe {
    let f_pv: *c_void = cast::transmute(f); 
    let gq_1D_integrand_caller_pv: *c_void = cast::transmute(gq_1D_integrand_caller);
    gauss_legendre(gq_order as c_int, gq_1D_integrand_caller_pv, f_pv, a, b)
  }
}




// This is the integrand callback function called directly by the space-adaptive C integration routine.
// Its job is to calculate the integrand value using f_ptr and set that value in *fval. 
#[inline(never)]
extern fn h_integrand_caller(ndim: c_uint, x: *R, 
                             f_ptr: *|x:&[R]| -> R, 
                             _: c_uint, fval: *mut R) -> c_int {
  unsafe {
    let f = ptr::read_ptr(f_ptr);
    *fval = vec::raw::buf_as_slice(x, ndim as uint, f);
  }
  0 as c_int
}

// This is the integrand callback function called directly by the C 2D gaussian quadrature integration routine.
#[inline(never)]
extern fn gq_2D_integrand_caller(x: R, y: R, f_ptr: *|x: R, y: R| -> R) -> R {
  let f: |x: R, y: R| -> R = unsafe { ptr::read_ptr(f_ptr) };
  f(x, y)
}

// This is the integrand callback function called directly by the C 1D gaussian quadrature integration routine.
#[inline(never)]
extern fn gq_1D_integrand_caller(x: R, f_ptr: *|R| -> R) -> R {
  let f: |R| -> R = unsafe { ptr::read_ptr(f_ptr) };
  f(x)
}



// The external C integration routines.
#[link_args = "lib/adaptive_quadrature.o lib/gaussian_quadrature.o"]
extern {

  fn hquadrature(fdim: c_uint,
                 f: *c_void, fdata: *c_void,
                 dim: c_uint, xmin: *R, xmax: *R, 
                 maxEval: size_t, reqAbsError: R, reqRelError: R, 
                 norm: u32,
                 val: *mut R, err: *mut R) -> c_int;

  fn gauss_legendre_2D_rect(order: c_int,
                            f: *c_void,    // double (*f)(double,double,void*)
                            data: *c_void, // passed as last param to f for each evaluation
                            a: R, b: R, c: R, d: R) -> R;
  
  fn gauss_legendre(order: c_int,
                    f: *c_void,    // double (*f)(double,void*)
                    data: *c_void, // passed as last param to f for each evaluation
                    a: R, b: R) -> R;
}


#[test]
fn test_h_quadrature() {
  let f1 = |_: &[f64]| 2.0;
  let f2 = |x: &[f64]| 2.0*x[0]*x[1];
  let min_bounds = ~[0.,0.];
  let max_bounds = ~[1.,1.];
  assert_eq!(space_adaptive_quadrature(&f1, min_bounds, max_bounds, 1e-5, 1e-5), 2.0)
  assert_eq!(space_adaptive_quadrature(&f2, min_bounds, max_bounds, 1e-5, 1e-5), 0.5)
  assert_eq!(space_adaptive_quadrature(&f1, min_bounds, max_bounds, 1e-5, 1e-5), 2.0)
}


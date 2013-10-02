use common::*;
use std::num::abs;

use polynomial::{Polynomial, PolyOwned};
use monomial::{Monomial, Mon2d, MaxMonDeg}; 
use mesh::{OShape};
use rectangle_mesh::{RectMesh, MeshCoord};
use weak_gradient::*;

// Tests below use the fact that if a polynomial p is in P_k(e) for some finite element e, where the
// weak gradient approximation space is [P_{k-1}(e)]^d, then wgrad(p) = grad(p) on the interior of e.
// Thus we can test weak gradient calculations on such polynomials and verify this equality. However,
// the weak gradient module only computes weak gradients of shape functions defined piecewise to be
// monomials on a single supporting face and 0 elsewhere. So to compute weak gradients of a polynomial
// across an entire finite element, we consider the polynomial as the sum of face-relative shape functions,
// and using the linearity of the weak gradient, we then compute the original polynomial's weak gradient
// as the sum of the weak gradients of these shape functions.

#[test]
fn test_wgrad_xy() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0f64, 0.],
                                              ~[3f64, 3.],
                                              ~[MeshCoord(3), MeshCoord(3)]);
  let wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(3), rmesh);

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
 
  let oshape_wgrads = wgrad_solver.wgrads_on_oshape([x*y],   // interior 
                                                    [&[one], // left side (ignored because xy is 0 here)
                                                     &[y],   // right
                                                     &[one], // bottom    (ignored ")
                                                     &[x]],  // top
                                                    OShape(0),
                                                    rmesh);

  let xy_on_int_wgrad = &oshape_wgrads.int_mon_wgrads[0];
  let y_on_right_side_wgrad = &oshape_wgrads.side_mon_wgrads[1][0]; 
  let x_on_top_side_wgrad = &oshape_wgrads.side_mon_wgrads[3][0];

  let wgrad = &WeakGrad::lcomb([(1., xy_on_int_wgrad), (1., y_on_right_side_wgrad), (1., x_on_top_side_wgrad)]);

  let wgrad_comp0 = PolyOwned::new(wgrad.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons.clone());
  let wgrad_comp1 = PolyOwned::new(wgrad.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons.clone());
 
  assert!(approx_equiv(&wgrad_comp0, &PolyOwned::new(~[1.],~[y]), 1e-9));
  assert!(approx_equiv(&wgrad_comp1, &PolyOwned::new(~[1.],~[x]), 1e-9));
}

#[test]
fn test_wgrad_x2y() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0f64, 0.],
                                              ~[3f64, 3.],
                                              ~[MeshCoord(3), MeshCoord(3)]);
  let wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(3), rmesh);

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
 
  let oshape_wgrads = wgrad_solver.wgrads_on_oshape([x*x*y],   // interior 
                                                    [&[one],   // left side (ignored because xy is 0 here)
                                                     &[y],     // right
                                                     &[one],   // bottom    (ignored ")
                                                     &[x*x]],  // top
                                                    OShape(0),
                                                    rmesh);

  let x2y_on_int_wgrad = &oshape_wgrads.int_mon_wgrads[0];
  let y_on_right_side_wgrad = &oshape_wgrads.side_mon_wgrads[1][0]; 
  let x2_on_top_side_wgrad = &oshape_wgrads.side_mon_wgrads[3][0];

  let wgrad = &WeakGrad::lcomb([(1., x2y_on_int_wgrad), (1., y_on_right_side_wgrad), (1., x2_on_top_side_wgrad)]);

  let wgrad_comp0 = PolyOwned::new(wgrad.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons.clone());
  let wgrad_comp1 = PolyOwned::new(wgrad.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons.clone());
 
  assert!(approx_equiv(&wgrad_comp0, &PolyOwned::new(~[2.],~[x*y]), 1e-9));
  assert!(approx_equiv(&wgrad_comp1, &PolyOwned::new(~[1.],~[x*x]), 1e-9));
}


fn approx_equiv<M:Monomial,P:Polynomial<M>>(p1: &P, p2: &P, tol: R) -> bool {
  let diff = p1 + p2.scaled(-1.);
  let diff_canon = diff.canonical_form();
  diff_canon.coefs.iter().all(|&c| abs(c) <= tol)
}


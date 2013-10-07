use weak_gradient::*;
use polynomial;
use polynomial::{Polynomial, PolyOwning, PolyBorrowingMons, approx_equiv};
use monomial::{Mon2d, MaxMonDeg}; 
use mesh::{OShape};
use rectangle_mesh::{RectMesh, MeshCoord};

use common::*;
use std::vec;

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
  let mut wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(3), rmesh);

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
 
  let (int_mon_wgrads, side_mon_wgrads) =
    wgrad_solver.wgrads_on_oshape([x*y],   // interior 
                                  [&[one], // left side (ignored because xy is 0 here)
                                  &[y],   // right
                                  &[one], // bottom    (ignored ")
                                  &[x]],  // top
                                  OShape(0),
                                  rmesh);

  let xy_on_int_wgrad = &int_mon_wgrads[0];
  let y_on_right_side_wgrad = &side_mon_wgrads[1][0]; 
  let x_on_top_side_wgrad = &side_mon_wgrads[3][0];

  let wgrad = &lcomb_wgrads([(1., xy_on_int_wgrad), (1., y_on_right_side_wgrad), (1., x_on_top_side_wgrad)]);

  let wgrad_comp0 = PolyOwning::new(wgrad.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons.clone());
  let wgrad_comp1 = PolyOwning::new(wgrad.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons.clone());
 
  assert!(approx_equiv(&wgrad_comp0, &PolyOwning::new(~[1.],~[y]), 1e-9));
  assert!(approx_equiv(&wgrad_comp1, &PolyOwning::new(~[1.],~[x]), 1e-9));
}

#[test]
fn test_wgrad_x2y() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0f64, 0.],
                                              ~[3f64, 3.],
                                              ~[MeshCoord(3), MeshCoord(3)]);
  let mut wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(3), rmesh);

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
 
  let (int_mon_wgrads, side_mon_wgrads) = 
    wgrad_solver.wgrads_on_oshape([x*x*y],   // interior 
                                  [&[one],   // left side (ignored because xy is 0 here)
                                  &[y],     // right
                                  &[one],   // bottom    (ignored ")
                                  &[x*x]],  // top
                                  OShape(0),
                                  rmesh);

  let x2y_on_int_wgrad = &int_mon_wgrads[0];
  let y_on_right_side_wgrad = &side_mon_wgrads[1][0]; 
  let x2_on_top_side_wgrad = &side_mon_wgrads[3][0];

  let wgrad = &lcomb_wgrads([(1., x2y_on_int_wgrad), (1., y_on_right_side_wgrad), (1., x2_on_top_side_wgrad)]);

  let wgrad_comp0 = PolyOwning::new(wgrad.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons.clone());
  let wgrad_comp1 = PolyOwning::new(wgrad.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons.clone());
 
  assert!(approx_equiv(&wgrad_comp0, &PolyOwning::new(~[2.],~[x*y]), 1e-9));
  assert!(approx_equiv(&wgrad_comp1, &PolyOwning::new(~[1.],~[x*x]), 1e-9));
}

#[test]
fn test_wgrad_dot() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0f64, 0.],
                                              ~[3f64, 3.],
                                              ~[MeshCoord(3), MeshCoord(3)]);
  let wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(2), rmesh);

  let wgrad1 = WeakGrad { comp_mon_coefs: ~[~[-2.5, 3., 0., 15., -6., -15.], ~[-2.5, 15., -15., 3., -6., 0.]] };
  let wgrad2 = WeakGrad { comp_mon_coefs: ~[~[2.5, -2., -0., -15., 6., 15.], ~[0., 0., 0., 0., 0., 0.]] }; 

  let mut wgrad_ops = wgrad_solver.weak_grad_ops();
  
  let fast_dot_prod = wgrad_ops.dot(&wgrad1, &wgrad2);
 
  // Compiler prevents second call of dot() which would have silently overwritten the existing wgrad's shared coefficients buffer. Beautiful!
  // let wgrad_x_top_dot_wgrad_y_right = wgrad_ops.dot(&wgrad1, &wgrad1);
  //   =>  error: cannot borrow `*wgrad_ops` as mutable more than once at a time
 
  let p1_0 = PolyBorrowingMons::new(wgrad1.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons);
  let p1_1 = PolyBorrowingMons::new(wgrad1.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons);
  let p2_0 = PolyBorrowingMons::new(wgrad2.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons);
  let p2_1 = PolyBorrowingMons::new(wgrad2.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons);
  let polys_dot_prod = PolyOwning::from_polys_lcomb([(1., &polynomial::mul(&p1_0, &p2_0)), (1., &polynomial::mul(&p1_1, &p2_1))]);

  println!("fast dot product:  {}", fast_dot_prod.to_str());
  println!("polys dot product: {}", polys_dot_prod.to_str());

  let canon_fast_dot_prod = fast_dot_prod.canonical_form();
  assert_eq!(&canon_fast_dot_prod.coefs, &polys_dot_prod.coefs);
  assert_eq!(&canon_fast_dot_prod.mons, &polys_dot_prod.mons);
}

fn lcomb_wgrads(terms: &[(R,&WeakGrad)]) -> WeakGrad {
  if terms.len() == 0 { fail!("lcomb_wgrads: At least one weak gradient is required.") }
  let (space_dims, num_comp_mons) = match terms[0] { (_, wgrad) => (wgrad.comp_mon_coefs.len(), wgrad.comp_mon_coefs[0].len()) };
  WeakGrad {
    comp_mon_coefs: 
      vec::from_fn(space_dims, |d| 
        vec::from_fn(num_comp_mons, |mon_num|
          terms.iter().fold(0 as R, |sum, &(c, wgrad)| sum + c * wgrad.comp_mon_coefs[d][mon_num])))
  }
}


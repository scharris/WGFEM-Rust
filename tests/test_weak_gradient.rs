use weak_gradient::*;
use polynomial;
use polynomial::{Polynomial, PolyOwning, PolyBorrowingMons, approx_equiv};
use monomial::{Mon2d, MaxMonDeg}; 
use mesh::{OShape};
use rectangle_mesh::{RectMesh, MeshCoord};
use dense_matrix::DenseMatrix;
use lapack;

use common::*;
use std::vec;

#[test]
fn test_do_lapack_init() {
  lapack::init(); // TODO: Do this somewhere else, where it's gauranteed to be run before other tests as part of each test setup.
}

#[test]
fn test_top_x_wgrad() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
  let mut wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(1), rmesh);

  let x = Mon2d { exps: [Deg(1), Deg(0)] };

  /* Component mon deg lim = 1
   * {q_i}_i = { (1,0), (y,0), (x,0), (0,1), (0,y), (0,x) }
   * v = x on top side, 0 elsewhere
   *
   * The system matrix contains integrals of the dot products of interior relative vector monomials over the interior
   * of the finite element, which in our case will be the unit square.
   * IP = inner products ({q_i}_i) = 
   *   int{(1,0).(1,0)}    int{(1,0).(y,0)}  int{(1,0).(x,0)} int{(1,0).(0,1)}, int{(1,0).(0,y)}  int{(1,0).(0,x)}
   *   int{(y,0).(1,0)}    int{(y,0).(y,0)}  int{(y,0).(x,0)} int{(y,0).(0,1)}, int{(y,0).(0,y)}  int{(y,0).(0,x)}
   *   int{(x,0).(1,0)}    int{(x,0).(y,0)}  int{(x,0).(x,0)} int{(x,0).(0,1)}, int{(x,0).(0,y)}  int{(x,0).(0,x)}
   *   int{(0,1).(1,0)}    int{(0,1).(y,0)}  int{(0,1).(x,0)} int{(0,1).(0,1)}, int{(0,1).(0,y)}  int{(0,1).(0,x)}
   *   int{(0,y).(1,0)}    int{(0,y).(y,0)}  int{(0,y).(x,0)} int{(0,y).(0,1)}, int{(0,y).(0,y)}  int{(0,y).(0,x)}
   *   int{(0,x).(1,0)}    int{(0,x).(y,0)}  int{(0,x).(x,0)} int{(0,x).(0,1)}, int{(0,x).(0,y)}  int{(0,x).(0,x)}
   *   =
   *            1                   1/2               1/2              0                 0                 0
   *           1/2                  1/3               1/4              0                 0                 0
   *           1/2                  1/4               1/3              0                 0                 0
   *            0                    0                 0               1                1/2               1/2
   *            0                    0                 0              1/2               1/3               1/4
   *            0                    0                 0              1/2               1/4               1/3
   * [WGRAD_DEF]
   *   (wgrad(v), q)_T = -(v_0, div q)_T + <v_b, q.n>_bnd(T), for all q in [P_r(T)]^d
   *
   * v_0 is 0 because our v is supported on the top side, so the first term of the RHS is 0 for each q.
   * [Beware:
   *   Here, as in the method itself, v is expressed as a *face-relative* monomial on its support face, while q is
   *   naturally expressed as an interior-relative vector monomial. This means that when v is side-supported, care
   *   must be taken to use a common coordinate system when for calculations involving q and v. For rectangle
   *   elements, this doesn't present any difficulty because the two coordinate systems will only differ in a
   *   single dimension (the side's perpendicular axis dimension), and monomials supported on a side face will
   *   be missing the problematic factor variable (otherwise they would be 0 and could not be part of a basis).
   *   For other element shapes conversion of at least one of q and v to a common coordinate system should be done
   *   to obtain correct results.]
   * RHS = 
   *   [q=(1,0)]: <x, (1,0).(0,1)>_top(T) = 0 
   *   [q=(y,0)]: <x, (y,0).(0,1)>_top(T) = 0  
   *   [q=(x,0)]: <x, (x,0).(0,1)>_top(T) = 0  
   *   [q=(0,1)]: <x, (0,1).(0,1)>_top(T) = 1/2
   *   [q=(0,y)]: <x, (0,y).(0,1)>_top(T) = 1/2 
   *   [q=(0,x)]: <x, (0,x).(0,1)>_top(T) = 1/3  
   *
   * Now we solve the system
   *   IP xi = RHS
   * This gives us the our solution as a linear combination of the q_i: 
   *   sol = sum_i xi_i q_i
   *   where
   *     xi = [0 0 0 -3/2 3 1]'
   */
  
  let (_, side_mon_wgrads) =
    wgrad_solver.wgrads_on_oshape([x],    // interior 
                                  [&[x],  // left side
                                   &[x],  // right
                                   &[x],  // bottom
                                   &[x]], // top
                                  OShape(0),
                                  rmesh);
  
  let x_on_top_side_wgrad = &side_mon_wgrads[3][0];

  assert_eq!(x_on_top_side_wgrad.comp_mon_coefs[0].as_slice(), &[0., 0., 0.]);
  assert_eq!(x_on_top_side_wgrad.comp_mon_coefs[1].as_slice(), &[-3./2., 3., 1.]);
}

#[test]
fn test_right_y_wgrad() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
  let mut wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(1), rmesh);

  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  /* System matrix is as above, only the right hand side of the system differs.
   *
   * [WGRAD_DEF]
   *   (wgrad(v), q)_T = -(v_0, div q)_T + <v_b, q.n>_bnd(T), for all q in [P_r(T)]^d
   *
   * RHS = 
   *   [q=(1,0)]: <y, (1,0).(1,0)>_top(T) = 1/2 
   *   [q=(y,0)]: <y, (y,0).(1,0)>_top(T) = 1/3 
   *   [q=(x,0)]: <y, (x,0).(1,0)>_top(T) = 1/2  
   *   [q=(0,1)]: <y, (0,1).(1,0)>_top(T) = 0
   *   [q=(0,y)]: <y, (0,y).(1,0)>_top(T) = 0
   *   [q=(0,x)]: <y, (0,x).(1,0)>_top(T) = 0
   *
   * Now we solve the system
   *   IP xi = RHS
   * This gives us the our solution as a linear combination of the q_i: 
   *   sol = sum_i xi_i q_i
   *   where
   *     xi = [-3/2 1 3 0 0 0]'
   */
  
  let (_, side_mon_wgrads) =
    wgrad_solver.wgrads_on_oshape([y],    // interior 
                                  [&[y],  // left side
                                   &[y],  // right
                                   &[y],  // bottom
                                   &[y]], // top
                                  OShape(0),
                                  rmesh);
  
  let y_on_right_side_wgrad = &side_mon_wgrads[1][0];

  assert_eq!(y_on_right_side_wgrad.comp_mon_coefs[0].as_slice(), &[-3./2., 1., 3.]);
  assert_eq!(y_on_right_side_wgrad.comp_mon_coefs[1].as_slice(), &[0., 0., 0.]);
}

#[test]
fn test_int_xy_wgrad() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
  let mut wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(1), rmesh);

  let xy = Mon2d { exps: [Deg(1), Deg(1)] };

  /* System matrix is as above, only the right hand side of the system differs.
   *
   * [WGRAD_DEF]
   *   (wgrad(v), q)_T = -(v_0, div q)_T + <v_b, q.n>_bnd(T), for all q in [P_r(T)]^d
   *
   * RHS = 
   *   [q=(1,0)]: -(xy, div (1,0))_T = 0
   *   [q=(y,0)]: -(xy, div (y,0))_T = 0
   *   [q=(x,0)]: -(xy, div (x,0))_T = -1/4
   *   [q=(0,1)]: -(xy, div (0,1))_T = 0
   *   [q=(0,y)]: -(xy, div (0,y))_T = -1/4
   *   [q=(0,x)]: -(xy, div (0,x))_T = 0
   *
   * Now we solve the system
   *   IP xi = RHS
   * This gives us the our solution as a linear combination of the q_i: 
   *   sol = sum_i xi_i q_i
   *   where
   *     xi = [3/2 0 -3 3/2 -3 0]'
   */
  
  let (int_mon_wgrads, _) =
    wgrad_solver.wgrads_on_oshape([xy],    // interior 
                                  [&[xy],  // left side
                                   &[xy],  // right
                                   &[xy],  // bottom
                                   &[xy]], // top
                                  OShape(0),
                                  rmesh);
  
  let xy_on_int_wgrad = &int_mon_wgrads[0];

  assert_eq!(xy_on_int_wgrad.comp_mon_coefs[0].as_slice(), &[3./2., 0., -3.]);
  assert_eq!(xy_on_int_wgrad.comp_mon_coefs[1].as_slice(), &[3./2., -3., 0.]);
}


// Tests below use the fact that if a polynomial p is in P_k(e) for some finite element e, where the
// weak gradient approximation space is [P_{k-1}(e)]^d, then wgrad(p) = grad(p) on the interior of e.
// Thus we can test weak gradient calculations on such polynomials and verify this equality. However,
// the weak gradient module only computes weak gradients of shape functions defined piecewise to be
// monomials on a single supporting face and 0 elsewhere. So to compute weak gradients of a polynomial
// across an entire finite element, we consider the polynomial as the sum of face-relative shape functions,
// and using the linearity of the weak gradient, we then compute the original polynomial's weak gradient
// as the sum of the weak gradients of these shape functions.

#[test]
fn test_wgrad_full_xy() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
  let mut wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(3), rmesh);

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
 
  let (int_mon_wgrads, side_mon_wgrads) =
    wgrad_solver.wgrads_on_oshape([x*y],    // interior 
                                  [&[one],  // left side (ignored because xy is 0 here)
                                   &[y],    // right
                                   &[one],  // bottom    (ignored ")
                                   &[x]],   // top
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
fn test_wgrad_full_x2y() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
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
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
  let wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(2), rmesh);

  let wgrad1 = WeakGrad { comp_mon_coefs: ~[~[-2.5, 3., 0., 15., -6., -15.], ~[-2.5, 15., -15., 3., -6., 0.]] };
  let wgrad2 = WeakGrad { comp_mon_coefs: ~[~[2.5, -2., -0., -15., 6., 15.], ~[0., 0., 0., 0., 0., 0.]] }; 

  let mut wgrad_ops = wgrad_solver.new_weak_grad_ops();
  
  let fast_dot_prod = wgrad_ops.dot(&wgrad1, &wgrad2);
 
  // Compiler prevents second call of dot() which would have silently overwritten the existing wgrad's shared coefficients
  // buffer (beautiful!).
  // let wgrad_x_top_dot_wgrad_y_right = wgrad_ops.dot(&wgrad1, &wgrad1);
  //   =>  error: cannot borrow `*wgrad_ops` as mutable more than once at a time

  // Effect a brute force dot product using polynomials for comparison.
  let p1_0 = PolyBorrowingMons::new(wgrad1.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons);
  let p1_1 = PolyBorrowingMons::new(wgrad1.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons);
  let p2_0 = PolyBorrowingMons::new(wgrad2.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons);
  let p2_1 = PolyBorrowingMons::new(wgrad2.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons);
  let polys_dot_prod = PolyOwning::from_polys_lcomb([(1., &polynomial::mul(&p1_0, &p2_0)),
                                                     (1., &polynomial::mul(&p1_1, &p2_1))]);

  let canon_fast_dot_prod = fast_dot_prod.canonical_form();
  assert_eq!(&canon_fast_dot_prod.coefs, &polys_dot_prod.coefs);
  assert_eq!(&canon_fast_dot_prod.mons,  &polys_dot_prod.mons);
}

#[test]
fn test_wgrad_mdot() {
  let rmesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0f64, 0.], ~[3f64, 3.], ~[MeshCoord(3), MeshCoord(3)]);
  let wgrad_solver: WeakGradSolver<Mon2d> = WeakGradSolver::new(MaxMonDeg(2), rmesh);

  // row swap matrix
  let m = DenseMatrix::from_fn(2,2, |r,c| { if r != c { 1. } else { 0. } });

  let wgrad1 = WeakGrad { comp_mon_coefs: ~[~[-1., 2., -3., 4., -5., 6.], ~[-2., 5., -3., 3., -6., 2.]] };
  let wgrad2 = WeakGrad { comp_mon_coefs: ~[~[ 5.,-4., -3., 2., 1., -2.], ~[ 1., 2.,  3., 4.,  5., 6.]] }; 

  let mut wgrad_ops = wgrad_solver.new_weak_grad_ops();
 
  let fast_dot_prod = wgrad_ops.mdot(&m, &wgrad1, &wgrad2);
 
  // Effect a brute force dot product using polynomials for comparison.
  let p1_0 = PolyBorrowingMons::new(wgrad1.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons);
  let p1_1 = PolyBorrowingMons::new(wgrad1.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons);
  let p2_0 = PolyBorrowingMons::new(wgrad2.comp_mon_coefs[0].clone(), wgrad_solver.wgrad_comp_mons);
  let p2_1 = PolyBorrowingMons::new(wgrad2.comp_mon_coefs[1].clone(), wgrad_solver.wgrad_comp_mons);
  let polys_dot_prod = PolyOwning::from_polys_lcomb([(1., &polynomial::mul(&p1_1, &p2_0)),   // p1's components swapped
                                                     (1., &polynomial::mul(&p1_0, &p2_1))]);

  let canon_fast_dot_prod = fast_dot_prod.canonical_form();
  assert_eq!(&canon_fast_dot_prod.coefs, &polys_dot_prod.coefs);
  assert_eq!(&canon_fast_dot_prod.mons,  &polys_dot_prod.mons);
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


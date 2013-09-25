use common::*;

use polynomial::{PolyWithBorrowedMons, poly};
use monomial::{Monomial, Mon2d}; 
use mesh::{OShape};
use rectangle_mesh::{RectMesh, MeshCoord, lesser_side_face_perp_to_axis, greater_side_face_perp_to_axis};
use weak_gradient::*;

fn test_wgrad_xy() {
  let rmesh: ~RectMesh<Mon2d> = RectMesh::new(~[0f64, 0.],
                                              ~[3f64, 3.],
                                              ~[MeshCoord(3), MeshCoord(3)]);
  let wgrad_solver: WeakGradientSolver<Mon2d> = WeakGradientSolver::new(Deg(3), rmesh);
  let left_side = lesser_side_face_perp_to_axis(Dim(0));
  let right_side = greater_side_face_perp_to_axis(Dim(0));
  let bottom_side = lesser_side_face_perp_to_axis(Dim(1));
  let top_side = greater_side_face_perp_to_axis(Dim(1));

  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };

  let wgrad = wgrad_solver.wgrad_int_mon(x*y, OShape(0), rmesh) +
              wgrad_solver.wgrad_side_mon(x, OShape(0), top_side, rmesh) +
              wgrad_solver.wgrad_side_mon(y, OShape(0), right_side, rmesh);
              // bottom and left are 0
  // TODO: check that wgrad is close to ~[poly([(1.,x)]),poly([(1.,y)])]
}


/*
zero = 0*one

# Test wgrad of monomials.

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2, rect_oshape, top_face, wgrad_solver) +
    wgrad(y, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y,1x^2])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y + 2.3, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2 + 2.3, rect_oshape, top_face, wgrad_solver) +
    wgrad(y + 2.3, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero + 2.3, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero + 2.3, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y,1x^2])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y^2 + -2., rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2 + -2., rect_oshape, top_face, wgrad_solver) +
    wgrad(y^2 + -2., rect_oshape, right_face, wgrad_solver) +
    wgrad(zero + -2., rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero + -2., rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y^2,2x^2*y])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y^2 + -2x*y + 1, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2 + -2x + 1, rect_oshape, top_face, wgrad_solver) +
    wgrad(y^2 + -2y + 1, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero + 1, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero + 1, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y^2 + -2y,2x^2*y + -2x])
)
*/


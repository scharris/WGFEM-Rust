use common::{R, Deg};
use monomial::{Mon2d, MaxMonDeg};
use mesh::{Mesh};
use rectangle_mesh::{RectMesh, MeshCoord};
use wg_basis::{WGBasis};
use vbf_laplace::VBFLaplace;
use wg_solver;
use wg_error_estimates::{err_L2_norm};
use lapack;

use std::num::{sin, cos};
use std::num::log2;
use std::iter::{range_step};

#[main]
fn main() {
  lapack::init();

  // u(x) = cos(x_0) + sin(x_1)
  // (grad u)(x) = (-sin(x_0), cos(x_1))
  // (div (grad u))(x) = -cos(x[0]) - sin(x[1])
  #[inline]
  fn u(x: &[R]) -> R { cos(x[0]) + sin(x[1]) }
  #[inline]
  fn f(x: &[R])-> R { cos(x[0]) + sin(x[1]) }
  #[inline]
  fn g(x: &[R])-> R { u(x) }

  //let poly_degs = range(2u8, 3);
  //let nums_side_divs = range_step(300u, 331, 10);
  let poly_degs = range(3u8, 4);
  let nums_side_divs = range_step(50u, 81, 10);

  let solve = |k: Deg, side_divs: uint| {
   
    let mesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[6.28,6.28], ~[MeshCoord(side_divs), MeshCoord(side_divs)]);
    let h = mesh.max_fe_diameter();
    println!("num fes: {}", mesh.num_fes());
    
    println("Constructing basis...");
    let basis = ~WGBasis::new(mesh, MaxMonDeg(*k), MaxMonDeg(*k-1));
    println!("num bels: {}", basis.num_els());
   
    let vbf = &VBFLaplace::new(None, basis);
 
    println("Solving system...");
    let wg_sol = wg_solver::solve(vbf, f, g);

    println("Computing L2 error...");
    let err = err_L2_norm(u, &wg_sol);
    println!("L2 error: {}", err);
    println("Done.");
    (h, err)
  };
  
  let deg_diam_L2errs: ~[(Deg,~[(R,R)])] = poly_degs.map(|k| {
    let diam_err_pairs = nums_side_divs.map(|side_divs| solve(Deg(k), side_divs)).collect();
    (Deg(k), diam_err_pairs)
  }).collect();

  for &(deg, ref diam_L2err_pairs) in deg_diam_L2errs.iter() {
    for i in range_step(0, diam_L2err_pairs.len(), 2) {
      if i + 1 < diam_L2err_pairs.len() {
        let ((h1,err1), (h2,err2)) = (diam_L2err_pairs[i], diam_L2err_pairs[i+1]);
        let slope = (log2(err2) - log2(err1))/(log2(h2) - log2(h1));
        println!("Deg {} mesh sizes {:f} and {:f}: convergence rate is {:f}", deg.to_str(), h1, h2, slope);
      }
    }
  }

}

/*
#[main]
fn main() {
  lapack::init();

  #[inline]
  fn u(x: &[R]) -> R { x[0]*(1.-x[0])*x[1]*(1.-x[1]) }
  #[inline]
  fn f(x: &[R])-> R { 2.*x[1]*(1.-x[1]) + 2.*x[0]*(1.-x[0]) }
  #[inline]
  fn g(x: &[R])-> R { 0 as R }

  let poly_degs = range(2u8, 3);
  let nums_side_divs = range_step(20u, 260, 10);

  let solve = |k: Deg, side_divs: uint| {
   
    let mesh: ~RectMesh<Mon2d> = ~RectMesh::new(~[0.,0.], ~[1.,1.], ~[MeshCoord(side_divs), MeshCoord(side_divs)]);
    let h = mesh.max_fe_diameter();
    println!("num fes: {}", mesh.num_fes());
    
    println("Constructing basis...");
    let basis = ~WGBasis::new(mesh, MaxMonDeg(*k), MaxMonDeg(*k-1));
    println!("num bels: {}", basis.num_els());
   
    let vbf = &VBFLaplace::new(None, basis);
 
    println("Solving system...");
    let wg_sol = wg_solver::solve(vbf, f, g);

    println("Computing L2 error...");
    let err = err_L2_norm(u, &wg_sol);
    println!("L2 error: {}", err);
    println("Done.");
    (h, err)
  };
  
  let deg_diam_L2errs: ~[(Deg,~[(R,R)])] = poly_degs.map(|k| {
    let diam_err_pairs = nums_side_divs.map(|side_divs| solve(Deg(k), side_divs)).collect();
    (Deg(k), diam_err_pairs)
  }).collect();

  for &(deg, ref diam_L2err_pairs) in deg_diam_L2errs.iter() {
    for i in range_step(0, diam_L2err_pairs.len(), 2) {
      if i + 1 < diam_L2err_pairs.len() {
        let ((h1,err1), (h2,err2)) = (diam_L2err_pairs[i], diam_L2err_pairs[i+1]);
        let slope = (log2(err2) - log2(err1))/(log2(h2) - log2(h1));
        println!("Deg {} mesh sizes {:f} and {:f}: convergence rate is {:f}", deg.to_str(), h1, h2, slope);
      }
    }
  }
}
*/


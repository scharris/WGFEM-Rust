use common::{R, sq};
use monomial::Monomial;
use mesh::{Mesh, FENum};
use wg_solution::WGSolution;

use std::num::sqrt;


pub fn err_L2_norm<Mon:Monomial,MeshT:Mesh<Mon>>(exact_sol: &fn(&[R])->R, approx: &WGSolution<Mon,MeshT>) -> R {
  let mesh = approx.basis().mesh();
  let sum_fe_sq_err_intgs = range(0, mesh.num_fes()).fold(0 as R, |sum_fe_sq_err_intgs, fe| {
    let sq_err_at_pt = |x: &[R], x_int_rel: &[R]| {
      sq(exact_sol(x) - approx.value_at_int_rel(FENum(fe), x_int_rel))
    };
    sum_fe_sq_err_intgs + mesh.intg_mixed_global_and_facerel_fn_on_fe_int(sq_err_at_pt, FENum(fe))
  });

  sqrt(sum_fe_sq_err_intgs)
}

/*
# Solution evaluation functions

# Return an interior-relative polynomial vector which represents the weak gradient of the passed
# solution on the indicated finite element interior.
function wg_sol_wgrad_on_interior(fe::FENum, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const wgrad_solver = basis.wgrad_solver
  const sol_coefs = wg_sol.basis_coefs
  const boundary_projs = wg_sol.boundary_projs
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)

  # The weak gradient of the solution is the sum of the weak gradients of the functions which equal
  # the solution on one face of our finite element and are zero elsewhere.
  total_wgrad = Poly.zero_poly_vec(Mesh.space_dim(mesh))

  # Add contributions for the terms of the solution's interior polynomial.
  for int_monn=monnum(1):basis.mons_per_fe_interior
    const beln = WGBasis.interior_mon_bel_num(fe, int_monn, basis)
    total_wgrad += sol_coefs[beln] * WGBasis.wgrad_interior_mon(int_monn, fe_oshape, basis)
  end

  # Add contributions for the solution's side polynomials, for both non-boundary and boundary sides.
  for sf=feface_one:Mesh.num_side_faces_for_shape(fe_oshape, mesh)
    if Mesh.is_boundary_side(fe, sf, mesh)
      const proj_coefs = boundary_projs[(fe, sf)]
      for side_monn=monnum(1):basis.mons_per_fe_side
        total_wgrad += proj_coefs[side_monn] * WGBasis.wgrad_side_mon(side_monn, fe_oshape, sf, basis)
      end
    else
      for side_monn=monnum(1):basis.mons_per_fe_side
        const beln = WGBasis.side_mon_bel_num(fe, sf, side_monn, basis)
        total_wgrad += sol_coefs[beln] * WGBasis.wgrad_side_mon(side_monn, fe_oshape, sf, basis)
      end
    end
  end

  total_wgrad
end


# Error Estimates


function err_H1_norm(exact_sol::Function,
                     grad_exact_sol::Array{Function,1},
                     wg_sol::WGSolution)
  const mesh = wg_sol.basis.mesh
  const d = Mesh.space_dim(mesh)

  # Work data for sum_sq_errs function below.
  const fe_rel_x = Array(R, d)
  const fe_origin = Array(R, d)

  integral_sq_errs = zeroR

  for fe=fenum(1):Mesh.num_fes(mesh)
    Mesh.fe_interior_origin!(fe, fe_origin, mesh)

    const approx = wg_sol_interior_poly(fe, wg_sol)
    const grad_approx = Poly.grad(approx)
    
    function sum_sq_errs(x::Vector{R})
      # set the fe-relative point coordinates
      for i=1:d
        fe_rel_x[i] = x[i] - fe_origin[i]
      end
      # compute the sum of the square first partial derivative errors
      sum_first_partial_sq_errs = zeroR
      for n=dim(1):d
        const partial_diff = grad_exact_sol[n](x) - Poly.value_at(grad_approx[n], fe_rel_x)
        sum_first_partial_sq_errs += partial_diff * partial_diff
      end
      const err = exact_sol(x) - Poly.value_at(approx, fe_rel_x)
      err * err + sum_first_partial_sq_errs
    end
    
    integral_sq_errs += Mesh.integral_global_x_face_rel_on_fe_face(sum_sq_errs,
                                                                   Mesh.one_mon(mesh),
                                                                   fe, Mesh.interior_face,
                                                                   mesh)
  end

  sqrt(integral_sq_errs)
end


function err_vs_proj_L2_norm(exact_sol::Function, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh

  integral_sq_err = zeroR
  for fe=fenum(1):Mesh.num_fes(mesh)
    const proj_poly = Proj.project_onto_fe_face_supported_approx_subspace_as_poly(exact_sol, fe, Mesh.interior_face, basis)
    const wg_sol_poly = wg_sol_interior_poly(fe, wg_sol)
    const err = proj_poly - wg_sol_poly
    const sq_err = Poly.canonical_form(err * err)
    const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
    integral_sq_err += Mesh.integral_face_rel_on_oshape_face(sq_err,
                                                             fe_oshape, Mesh.interior_face,
                                                             mesh)
  end

  sqrt(integral_sq_err)
end

function err_grad_vs_wgrad_L2_norm(exact_sol_grad::Function, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)
  const fe_origin = Array(R, d)

  integral_sq_err = zeroR
  for fe=fenum(1):Mesh.num_fes(mesh)
    Mesh.fe_interior_origin!(fe, fe_origin, mesh)
    const wgrad = wg_sol_wgrad_on_interior(fe, wg_sol)
    function sq_err(x::Vector{R})
      for i=1:d
        fe_rel_x[i] = x[i] - fe_origin[i]
      end
      const diff = exact_sol_grad(x) - Poly.value_at(wgrad, fe_rel_x)
      dot(diff, diff)
    end
    integral_sq_err += Mesh.integral_global_x_face_rel_on_fe_face(sq_err,
                                                                  one,
                                                                  fe, Mesh.interior_face,
                                                                  mesh)
  end
  sqrt(integral_sq_err)
end

# Compute semi-norm ||| Q_h u - u_h |||_vbf = sqrt(vbf(Q u - u_h, Q u - u_h)),
# where Q_h is piecwise L2 projection onto finite element faces.
function err_vs_proj_vbf_seminorm(exact_sol::Function, wg_sol::WGSolution, vbf::AbstractVariationalBilinearForm)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const vbf_symm = VBF.is_symmetric(vbf)

  # work arrays
  const side_diffs = Array(Vector{R}, Mesh.max_num_shape_sides(basis.mesh))
  const is_nb_side = Array(Bool, Mesh.max_num_shape_sides(basis.mesh))

  # We must evaluate
  #   sum_{fe} vbf(diff(fe), diff(fe)),
  # where
  #   diff(fe) = diff_int(fe) + sum_{s in nb-sides(fe)} diff_s,
  # and diff_f for a face f is the L2 projection onto f of the exact solution minus the wg
  # solution on the face, and 0 elsewhere.
  # By the bilinearity of the vbf, we need to sum the contributions from each finite element
  # of the vbf applied to the differences on each interior or non-boundary side against that
  # on each interior or non-boundary side.

  sum = zeroR
  for fe=fenum(1):Mesh.num_fes(mesh)
    const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
    const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, mesh)
    # Find non-boundary sides
    for sf=feface_one:num_sides
      is_nb_side[sf] = !Mesh.is_boundary_side(fe, sf, mesh)
    end

    # Compute exact solution projection minus wg solution for all faces.
    const int_diff = Proj.project_onto_fe_face_supported_approx_subspace(exact_sol, fe, Mesh.interior_face, basis) -
                     wg_sol_interior_coefs(fe, wg_sol)
    for sf=feface_one:num_sides if is_nb_side[sf]
      diff = Proj.project_onto_fe_face_supported_approx_subspace(exact_sol, fe, sf, basis)
      diff -= wg_sol_side_coefs(fe, sf, wg_sol)
      side_diffs[sf] = diff
    end end

    # interior vs interior
    sum += VBF.poly_on_face_vs_poly_on_face(fe_oshape,
                                            int_diff, Mesh.interior_face,
                                            int_diff, Mesh.interior_face,
                                            basis,
                                            vbf)

    for sf=feface_one:num_sides if is_nb_side[sf]
      const sf_diff = side_diffs[sf]

      # Add contributions from pairings of a side and interior.

      # interior vs side
      const int_vs_side = VBF.poly_on_face_vs_poly_on_face(fe_oshape,
                                                           int_diff, Mesh.interior_face,
                                                           sf_diff, sf,
                                                           basis,
                                                           vbf)
      sum += int_vs_side

      # side vs interior
      sum += if vbf_symm  int_vs_side
             else VBF.poly_on_face_vs_poly_on_face(fe_oshape,
                                                   sf_diff, sf,
                                                   int_diff, Mesh.interior_face,
                                                   basis,
                                                   vbf) end

      # Add contributions for pairings of sides.

      # side vs side
      if vbf_symm
        # contributions from pairings of unequal sides (symmetric case)
        for sf2=fefacenum(sf+1):num_sides if is_nb_side[sf2]
          sum += 2 * VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, side_diffs[sf2], sf2, basis, vbf)
        end end
        # contributions from side self-pairings (symmetric case)
        sum += VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, sf_diff, sf, basis, vbf)
      else # vbf not symmetric
        # contributions for all side vs side pairings (asymmetric case)
        for sf2=feface_one:num_sides if is_nb_side[sf2]
          sum += VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, side_diffs[sf2], sf2, basis, vbf)
        end end
      end
    end end # non-boundary sides
  end
  sqrt(sum)
end

end
*/

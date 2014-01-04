#[crate_id = "wgfem#0.8.1"];

#[license = "MIT"];
#[crate_type = "exe"];
#[feature(globs)];
#[feature(macro_rules)];
#[feature(link_args)];

extern mod extra;

pub mod common;
pub mod monomial;
pub mod polynomial;
pub mod vector_monomial;
mod quadrature;
#[path="linear_algebra.rs"] mod la;
pub mod dense_matrix;
pub mod sparse_matrix;
pub mod storage_by_ints;
pub mod mesh;
pub mod rectangle_mesh;
pub mod triangle_mesh;
pub mod triangle_mesh_builder;
pub mod weak_gradient;
pub mod wg_basis;
pub mod projection;
pub mod variational_bilinear_form;
pub mod vbf_laplace;
pub mod wg_solution;
pub mod wg_solver;
pub mod wg_error_estimates;
pub mod main;

#[cfg(test)]
mod tests {
  // no tests for common
  mod test_monomial;
  mod test_polynomial;
  mod test_vector_monomial;
  mod test_dense_matrix;
  mod test_sparse_matrix;
  mod test_storage_by_ints;
  mod test_la;
  // no tests for quadrature
  // no tests for abstract mesh
  mod test_rectangle_mesh;
  mod test_weak_gradient;
  mod test_wg_basis;
  mod test_projection;
  mod test_variational_bilinear_form;
  mod test_vbf_laplace;
}


#[link(name = "wgfem",
  vers = "0.1",
  url = "https://github.com/scharris/WGFEM-Rust")];

#[license = "MIT"];
#[crate_type = "lib"];
#[feature(globs)];
#[feature(macro_rules)];

extern mod extra;

pub mod common;
pub mod monomial;
pub mod polynomial;
pub mod vector_monomial;
pub mod mesh;
pub mod rectangle_mesh;
pub mod dense_matrix;
pub mod weak_gradient;
pub mod wg_basis;
pub mod projection;
pub mod sparse_matrix;
pub mod tensor;
pub mod variational_bilinear_form;
pub mod lapack;
mod quadrature;

#[cfg(test)]
mod tests {
  mod test_monomial;
  mod test_polynomial;
  mod test_vector_monomial;
  // no tests for abstract mesh
  mod test_rectangle_mesh;
  mod test_dense_matrix;
  mod test_sparse_matrix;
  mod test_tensor;
  mod test_weak_gradient;
  mod test_wg_basis;
  mod test_projection;
  mod test_variational_bilinear_form;
}


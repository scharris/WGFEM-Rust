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
#[cfg(test)]
mod test_monomial;

pub mod polynomial;
#[cfg(test)]
mod test_polynomial;

pub mod vector_monomial;
#[cfg(test)]
mod test_vector_monomial;

pub mod mesh;

pub mod rectangle_mesh;
#[cfg(test)]
mod test_rectangle_mesh;

pub mod dense_matrix;
#[cfg(test)]
mod test_dense_matrix;

pub mod weak_gradient;
#[cfg(test)]
mod test_weak_gradient;

pub mod wg_basis;
#[cfg(test)]
mod test_wg_basis;

pub mod lapack;
mod quadrature;


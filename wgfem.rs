#[link(name = "wgfem",
  vers = "0.1",
  url = "https://github.com/scharris/WGFEM-Rust")];

#[license = "MIT"];
#[crate_type = "lib"];

extern mod extra;

pub mod common;
pub mod monomial;
pub mod polynomial;
pub mod vector_monomial;
pub mod mesh;
pub mod rect_mesh;
pub mod dense_matrix;
pub mod weak_gradient;
mod quadrature;
#[cfg(test)]
mod test_monomial;
#[cfg(test)]
mod test_vector_monomial;
#[cfg(test)]
mod test_polynomial;
#[cfg(test)]
mod test_rect_mesh;
#[cfg(test)]
mod test_weak_gradient;


use common::*;
use vector_monomial::VectorMonomial;
use monomial::Monomial;
use mesh::Mesh;
use dense_matrix::DenseMatrix;

/*
 * For a weak function v on a finite element T, the weak gradient of degree r
 * of v on T is defined to be the polynomial vector function wgrad(v) in
 * [P_r(T)]^d, such that:
 *
 * [WGRAD_DEF]
 *   (wgrad(v), q)_T = -(v_0, div q)_T + <v_b, q.n>_bnd(T), for all q in [P_r(T)]^d
 *
 * By linearity of both sides in the right sides of the inner products, the above
 * holds iff the same equation holds with the q's further restricted to be
 * vector monomials on T (which are a monomial in one range component and 0 for
 * all other components), which functions form a basis {b_i}_i for [P_r(T)]^d.
 *
 * Letting q = b_i, and writing wgrad(v) = sum_j eta_j b_j, we obtain a linear
 * system from WGRAD_DEF which we can solve for the unkowns eta_j, with
 * (b_i, b_j)_T being the matrix elements of the system, and the right hand side
 * of WGRAD_DEF defining the right hand side column vector for the system.
 *
 * We will only actually need weak gradients of monomials or polynomials
 * supported on a single face (interior or side) of a finite element. Thus v
 * will be expressed below as a monomial or polynomial paired with the face of T
 * on which it has the monomial or polynomial value.
 */

struct WeakGradientSolver<M> {

  basis_vmons: ~[VectorMonomial<M>],

  ips_basis_vs_basis_by_oshape: ~[DenseMatrix]

  // NOTE: basis_divs
  // Probably don't need to store these divergences of the basis vmons, since they are so easily computed.
  // If the vmon is used in the same calculation, probably best to just compute the div instead of doing a fetch
  // of a precomputed value.
  // basis_divs: ~[(R,M)]

  // NOTE: mesh
  // Not necessary to store this in the structure, callers should just pass it in when needed.
}

impl <M:Monomial> WeakGradientSolver<M> {

  pub fn new<MESH:Mesh<M>>(wgrad_mons_max_deg: Deg, mesh: MESH) -> WeakGradientSolver<M> {
    let vmons: ~[VectorMonomial<M>] = VectorMonomial::vector_mons_of_deg_le(wgrad_mons_max_deg);
    let ips = WeakGradientSolver::make_vmon_ips_by_oshape(vmons, mesh);
    WeakGradientSolver { 
      basis_vmons: vmons,
      ips_basis_vs_basis_by_oshape: ips
    }
  }

  fn make_vmon_ips_by_oshape<MESH:Mesh<M>>(basis_vmons: &[VectorMonomial<M>], mesh: MESH) -> ~[DenseMatrix] {
    // TODO
    ~[]
  }
}


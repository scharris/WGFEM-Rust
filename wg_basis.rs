use common::*;
use monomial::{Monomial, DegLim, MaxMonDeg, MaxMonFactorDeg, domain_space_dims};
use mesh::{Mesh, OShape, SideFace};
use dense_matrix::DenseMatrix;
use weak_gradient::{WeakGradSolver, WeakGrad};

use std::vec;

// TODO: Discuss role of the basis type here and general strategy of enumeration.

// Side Monomials and Side Dependent Dimensions
// --------------------------------------------
// Side supported shape function sequences are stored by their declared side dependent dimension, which determines
// which monomials satisfying the degree limit can be included to avoid having a linearly dependent set. A side has
// dependent dimension of r if coordinate r can be expressed as a function (affine) of the other coordinates on the
// side. The dependent dimension for a side is determined by the mesh, and usually chosen from among multiple 
// acceptable alternatives. A side declared to have dependent dimension r by the mesh will host shape functions
// whose local monomials on the side are constant (have exponent 0) in the r^th coordinate factor.


#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct BasisElNum(uint);

#[deriving(Eq,TotalEq,Ord,TotalOrd,Clone)]
pub struct FaceMonNum(uint);


// A type representing a basis for Weak Galerkin approximating polynomials on an arbitrary mesh.
struct WgBasis<Mon,Mesh> {
  
  // The mesh over which the Weak Galerkin basis is formed.
  mesh: ~Mesh,

  // Degree limits of the approximating polynomials.
  int_polys_deg_lim: DegLim,
  side_polys_deg_lim: DegLim,

  // Monomial sequence used to define basis shape functions on any one finite element interior.
  int_mons: ~[Mon],

  // This vector stores the sequences of monomials used to define basis shape functions supported on 
  // any single finite element side, arranged by side dependent dimension. Thus in the vector at
  // position r are the sequence of side-local monomials defining the shape functions in the basis
  // which are supported on any given side which the mesh declares to have dependent dimension r.
  side_mons_by_dep_dim: ~[~[Mon]], // side monomial sequences by side dependent dimension

  mons_per_fe_int: uint,
  mons_per_fe_side: uint,

  // Significant counts in the basis element enumeration.
  total_els: uint,
  num_int_els: uint, 
  
  // First side supported basis element. All preceeding basis elements are interior supported.
  first_nb_side_beln: BasisElNum,

  // Weak gradients generator.
  weak_grad_solver: ~WeakGradSolver<Mon>,

  // Pre-calculated weak gradients of basis elements supported on reference oriented shapes.
  int_mon_wgrads: ~[~[WeakGrad]],     // by fe oshape, then interior monomial number
  side_mon_wgrads: ~[~[~[WeakGrad]]], // by fe oshape, then side face, then side monomial number

  // Pre-calculated L2 inner products between basis elements supported on the same faces of reference oriented shapes.
  // Only the upper triangle part of each matrix should be used, the contents of the lower parts are undefined.
  // These are provided for clients performing L2 projections onto fe interiors and sides.
  ips_int_mons: ~[DenseMatrix],     // by fe oriented shape, then (int mon num, int mon num)
  ips_side_mons: ~[~[DenseMatrix]], // by fe oriented shape, then side face, then (side mon #, side mon #)
}


impl <Mon:Monomial, MeshT:Mesh<Mon>> WgBasis<Mon,MeshT> {

  pub fn new(mesh: ~MeshT, int_polys_deg_lim: DegLim, side_polys_deg_lim: DegLim) -> ~WgBasis<Mon,MeshT> {
    
    let space_dims = domain_space_dims::<Mon>();
    
    let int_mons = Monomial::mons_with_deg_lim_asc(int_polys_deg_lim);
    
    let side_mons_by_dep_dim: ~[~[Mon]] = { 
      let mons_for_deg_lim: ~[Mon] = Monomial::mons_with_deg_lim_asc(side_polys_deg_lim);
      vec::from_fn(space_dims, |r|
        mons_for_deg_lim.iter().filter(|mon| mon.exp(Dim(r)) == Deg(0)).map(|m|m.clone()).collect()
      )
    };
    
    let mons_per_fe_int = int_mons.len();
    let mons_per_fe_side = side_mons_by_dep_dim[0].len();

    let num_int_els = mesh.num_fes() * int_mons.len(); 
    let total_els = num_int_els + mesh.num_nb_sides() * mons_per_fe_side;
    let first_nb_side_beln = BasisElNum(num_int_els + 1);

    let mut wgrad_solver = {
      let k = match int_polys_deg_lim { MaxMonDeg(k) => k, MaxMonFactorDeg(k) => k };
      ~WeakGradSolver::new(MaxMonDeg(k-1), mesh)
    };

    let (int_mon_wgrads, side_mon_wgrads) = compute_wgrads(wgrad_solver, int_mons, side_mons_by_dep_dim, mesh);
    
    let ips_int_mons = make_int_mon_ips(int_mons, mesh);
    let ips_side_mons = make_side_mon_ips(side_mons_by_dep_dim, mesh);

    ~WgBasis {
      mesh: mesh,
      int_polys_deg_lim: int_polys_deg_lim,
      side_polys_deg_lim: side_polys_deg_lim,
      int_mons: int_mons,
      side_mons_by_dep_dim: side_mons_by_dep_dim,
      mons_per_fe_int: mons_per_fe_int,
      mons_per_fe_side: mons_per_fe_side,
      total_els: total_els,
      num_int_els: num_int_els,
      first_nb_side_beln: first_nb_side_beln,
      weak_grad_solver: wgrad_solver,
      int_mon_wgrads: int_mon_wgrads,
      side_mon_wgrads: side_mon_wgrads,
      ips_int_mons: ips_int_mons,
      ips_side_mons: ips_side_mons,
    }
  }

}

fn compute_wgrads<Mon:Monomial,MeshT:Mesh<Mon>>(wgrad_solver: &mut WeakGradSolver<Mon>,
                                                int_mons: &[Mon],
                                                side_mons_by_dep_dim: &[~[Mon]],
                                                mesh: &MeshT) -> (~[~[WeakGrad]], ~[~[~[WeakGrad]]]) {
  let mut int_mon_wgrads_by_oshape = vec::with_capacity(mesh.num_oriented_element_shapes());
  let mut side_mon_wgrads_by_oshape = vec::with_capacity(mesh.num_oriented_element_shapes());

  for os in range(0, mesh.num_oriented_element_shapes()) {
    let os = OShape(os);
    let side_mons_by_side = vec::from_fn(mesh.num_side_faces_for_shape(os), |sf| {
      let sf_dep_dim = mesh.dependent_dim_for_oshape_side(os, SideFace(sf));
      side_mons_by_dep_dim[*sf_dep_dim].as_slice()
    });
    
    let (int_mon_wgrads, side_mon_wgrads) = wgrad_solver.wgrads_on_oshape(int_mons, side_mons_by_side, os, mesh);

    int_mon_wgrads_by_oshape.push(int_mon_wgrads);
    side_mon_wgrads_by_oshape.push(side_mon_wgrads);
  } 

  (int_mon_wgrads_by_oshape, side_mon_wgrads_by_oshape)
}

fn make_int_mon_ips<Mon:Monomial,MeshT:Mesh<Mon>>(int_mons: &[Mon], mesh: &MeshT) -> ~[DenseMatrix] {
  vec::from_fn(mesh.num_oriented_element_shapes(), |os| {
    DenseMatrix::with_upper_triangle_from_fn(int_mons.len(), int_mons.len(), |i,j| {
      mesh.intg_facerel_mon_on_oshape_int(int_mons[i] * int_mons[j], OShape(os))
    })
  })
}

fn make_side_mon_ips<Mon:Monomial,MeshT:Mesh<Mon>>(side_mons_by_dep_dim: &[~[Mon]], mesh: &MeshT) -> ~[~[DenseMatrix]] {
  vec::from_fn(mesh.num_oriented_element_shapes(), |os| {
    vec::from_fn(mesh.num_side_faces_for_shape(OShape(os)), |sf| {
      let sf_dep_dim = mesh.dependent_dim_for_oshape_side(OShape(os), SideFace(sf));
      let side_mons = side_mons_by_dep_dim[*sf_dep_dim].as_slice();
      DenseMatrix::with_upper_triangle_from_fn(side_mons.len(), side_mons.len(), |i,j| {
        mesh.intg_facerel_mon_on_oshape_side(side_mons[i] * side_mons[j], OShape(os), SideFace(sf))
      })
    })
  })
}


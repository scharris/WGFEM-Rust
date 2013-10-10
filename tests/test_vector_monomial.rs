use common::*;
use monomial::{Monomial, Mon2d, Mon3d, Mon4d, MaxMonDeg};
use vector_monomial::*;

#[test]
fn test_construction() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let _: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x);
}

#[test]
#[should_fail]
fn test_improper_construction() {
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let _: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(3), x);
}

#[test]
fn test_divergence() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  
  let x_dim0_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x);
  assert_eq!(x_dim0_vmon.divergence_coef_and_mon(), (1.,one));
  
  let y_dim0_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), y);
  assert_eq!(y_dim0_vmon.divergence_coef_and_mon(), (0.,y));
  
  let y_dim1_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), y);
  assert_eq!(y_dim1_vmon.divergence_coef_and_mon(), (1.,one));
  
  let xy_dim0_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x*y);
  assert_eq!(xy_dim0_vmon.divergence_coef_and_mon(), (1.,y));

  let x2y3_dim1_vmon: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x*x*y*y*y);
  assert_eq!(x2y3_dim1_vmon.divergence_coef_and_mon(), (3., x*x*y*y));
}

#[test]
fn test_ordered_by_comp_mon_of_deg_le_2d() {
  let one = Mon2d { exps: [Deg(0), Deg(0)] };
  let x = Mon2d { exps: [Deg(1), Deg(0)] };
  let y = Mon2d { exps: [Deg(0), Deg(1)] };
  let one_dim0: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), one);
  let y_dim0: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), y);
  let y2_dim0: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), y*y);
  let x_dim0: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x);
  let xy_dim0: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x*y);
  let x2_dim0: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(0), x*x);
  let one_dim1: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), one);
  let y_dim1: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), y);
  let y2_dim1: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), y*y);
  let x_dim1: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x);
  let xy_dim1: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x*y);
  let x2_dim1: VectorMonomial<Mon2d> = VectorMonomial::new(Dim(1), x*x);

  let mons_deg0 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(0));
  let vmons_deg0: ~[VectorMonomial<Mon2d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg0);
  assert_eq!(&vmons_deg0, &~[one_dim0, one_dim1]);

  let mons_deg1 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(1));
  let vmons_deg1: ~[VectorMonomial<Mon2d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg1);
  assert_eq!(&vmons_deg1,
             &~[one_dim0, y_dim0, x_dim0, 
                one_dim1, y_dim1, x_dim1]);
  
  let mons_deg2 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(2));
  let vmons_deg2: ~[VectorMonomial<Mon2d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg2);
  assert_eq!(&vmons_deg2,
             &~[one_dim0, y_dim0, y2_dim0, x_dim0, xy_dim0, x2_dim0,
                one_dim1, y_dim1, y2_dim1, x_dim1, xy_dim1, x2_dim1]);
}

#[test]
fn test_ordered_by_comp_mon_of_deg_le_3d() {
  let one = Mon3d { exps: [Deg(0), Deg(0), Deg(0)] };
  let x = Mon3d { exps: [Deg(1), Deg(0), Deg(0)] };
  let y = Mon3d { exps: [Deg(0), Deg(1), Deg(0)] };
  let z = Mon3d { exps: [Deg(0), Deg(0), Deg(1)] };
  let one_dim0: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(0), one);
  let z_dim0: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(0), z);
  let y_dim0: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(0), y);
  let x_dim0: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(0), x);
  let one_dim1: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(1), one);
  let z_dim1: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(1), z);
  let y_dim1: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(1), y);
  let x_dim1: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(1), x);
  let one_dim2: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(2), one);
  let z_dim2: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(2), z);
  let y_dim2: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(2), y);
  let x_dim2: VectorMonomial<Mon3d> = VectorMonomial::new(Dim(2), x);
  
  let mons_deg0 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(0));
  let vmons_deg0: ~[VectorMonomial<Mon3d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg0);
  assert_eq!(&vmons_deg0, &~[one_dim0, one_dim1, one_dim2]);

  let mons_deg1 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(1));
  let vmons_deg1: ~[VectorMonomial<Mon3d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg1);
  assert_eq!(&vmons_deg1,
             &~[one_dim0, z_dim0, y_dim0, x_dim0, 
                one_dim1, z_dim1, y_dim1, x_dim1,
                one_dim2, z_dim2, y_dim2, x_dim2]);
}

#[test]
fn test_ordered_by_comp_mon_of_deg_le_4d() {
  let one = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(0)] };
  let x = Mon4d { exps: [Deg(1), Deg(0), Deg(0), Deg(0)] };
  let y = Mon4d { exps: [Deg(0), Deg(1), Deg(0), Deg(0)] };
  let z = Mon4d { exps: [Deg(0), Deg(0), Deg(1), Deg(0)] };
  let t = Mon4d { exps: [Deg(0), Deg(0), Deg(0), Deg(1)] };
  let one_dim0: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(0), one);
  let t_dim0: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(0), t);
  let z_dim0: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(0), z);
  let y_dim0: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(0), y);
  let x_dim0: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(0), x);
  let one_dim1: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(1), one);
  let t_dim1: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(1), t);
  let z_dim1: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(1), z);
  let y_dim1: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(1), y);
  let x_dim1: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(1), x);
  let one_dim2: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(2), one);
  let t_dim2: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(2), t);
  let z_dim2: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(2), z);
  let y_dim2: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(2), y);
  let x_dim2: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(2), x);
  let one_dim3: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(3), one);
  let t_dim3: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(3), t);
  let z_dim3: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(3), z);
  let y_dim3: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(3), y);
  let x_dim3: VectorMonomial<Mon4d> = VectorMonomial::new(Dim(3), x);
  
  let mons_deg0 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(0));
  let vmons_deg0: ~[VectorMonomial<Mon4d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg0);
  assert_eq!(&vmons_deg0, &~[one_dim0, one_dim1, one_dim2, one_dim3]);

  let mons_deg1 = Monomial::mons_with_deg_lim_asc(MaxMonDeg(1));
  let vmons_deg1: ~[VectorMonomial<Mon4d>] = VectorMonomial::with_comp_mons_ordered_by_comp_and_mon(mons_deg1);
  assert_eq!(&vmons_deg1,
             &~[one_dim0, t_dim0, z_dim0, y_dim0, x_dim0, 
                one_dim1, t_dim1, z_dim1, y_dim1, x_dim1,
                one_dim2, t_dim2, z_dim2, y_dim2, x_dim2,
                one_dim3, t_dim3, z_dim3, y_dim3, x_dim3]);
}


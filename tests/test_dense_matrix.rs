use dense_matrix::*;
use common::R;
use lapack;

use std::vec;

#[test]
fn test_do_lapack_init() {
  lapack::init(); // TODO: Do this somewhere else, where it's gauranteed to be run before other tests as part of each test setup.
}

#[test]
fn test_constr_from_fn() {
  let m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  
  let data_as_vec = unsafe { vec::from_buf(m.col_maj_data_ptr(), 9) };
  assert_eq!(data_as_vec, ~[0.,1000.,2000.,1.,1001.,2001.,2.,1002.,2002.]);
}

#[test]
fn test_constr_from_elem_with_cols_capacity() {
  let m = DenseMatrix::from_elem_with_cols_capacity(3,3, 10., 5);
  let data_as_vec = unsafe { vec::from_buf(m.col_maj_data_ptr(), 15) };
  assert_eq!(data_as_vec, ~[10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.]);
  assert_eq!(m.get(1,1), 10.);
  assert_eq!(m.get(2,2), 10.);
}

#[test]
fn test_constr_of_size_with_cols_capacity() {
  let mut m = DenseMatrix::of_size_with_cols_capacity(3,3,5);
  m.set(1,1, 11.);
  m.set(2,2, 22.);
  assert_eq!(m.get(1,1), 11.);
  assert_eq!(m.get(2,2), 22.);
}

#[test]
fn test_constr_from_rows() {
  let m = DenseMatrix::from_rows(3,3,
    [~[   0.,    1.,    2.],
     ~[1000., 1001., 1002.],
     ~[2000., 2001., 2002.]]);

  let data_as_vec = unsafe { vec::from_buf(m.col_maj_data_ptr(), 9) };
  assert_eq!(data_as_vec, ~[0.,1000.,2000.,1.,1001.,2001.,2.,1002.,2002.]);
}


#[test]
#[should_fail]
fn test_bad_col_access_under_capacity1() {
  let m = DenseMatrix::from_elem_with_cols_capacity(3,3, 10., 5);
  m.get(3,3);
}

#[test]
#[should_fail]
fn test_bad_col_access_under_capacity2() {
  let m = DenseMatrix::of_size_with_cols_capacity(3,3, 5);
  m.get(3,3);
}

#[test]
fn test_constr_upper_triangle_from_fn() {
  let m = DenseMatrix::upper_triangle_from_fn(3, |r,c| {
    1000. * (r as R) + (c as R)
  });
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
fn test_constr_lower_triangle_from_fn() {
  let m = DenseMatrix::lower_triangle_from_fn(3, |r,c| {
    1000. * (r as R) + (c as R)
  });
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(1,0), 1000.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(2,0), 2000.);
  assert_eq!(m.get(2,1), 2001.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
fn test_get() {
  let m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,0), 1000.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,0), 2000.);
  assert_eq!(m.get(2,1), 2001.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
fn test_set() {
  let mut m = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  assert_eq!(m.get(2,1), 2001.);
  m.set(2,1, 2001.5);
  assert_eq!(m.get(2,1), 2001.5);
}

#[test]
fn test_copy_into() {
  let m_src = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  let m = &mut DenseMatrix::from_elem(3,3, 10.);
  m_src.copy_into(m);
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,0), 1000.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,0), 2000.);
  assert_eq!(m.get(2,1), 2001.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
#[should_fail]
fn test_bad_copy_into() {
  let m_src = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  let m = &mut DenseMatrix::from_elem(3,2, 10.);
  m_src.copy_into(m);
}

#[test]
fn test_copy_upper_triangle_into() {
  let m_src = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  let m = &mut DenseMatrix::from_elem(3,3, 10.);
  m_src.copy_upper_triangle_into(m);
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
fn test_fill_upper_triangle_from() {
  let m_src = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  let m = &mut DenseMatrix::from_elem(3,3, 10.);
  m.fill_upper_triangle_from(&m_src);
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
fn test_fill_from_fn() {
  let mut m = DenseMatrix::from_elem(3,3, 0 as R);
  m.fill_from_fn(|r,c| 1000. * r as R + c as R);
  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(0,2), 2.);
  assert_eq!(m.get(1,0), 1000.);
  assert_eq!(m.get(1,1), 1001.);
  assert_eq!(m.get(1,2), 1002.);
  assert_eq!(m.get(2,0), 2000.);
  assert_eq!(m.get(2,1), 2001.);
  assert_eq!(m.get(2,2), 2002.);
}

#[test]
#[should_fail]
fn test_bad_copy_upper_triangle_into() {
  let m_src = DenseMatrix::from_fn(3,3, |r,c| {
    1000. * r as R + c as R
  });
  let m = &mut DenseMatrix::from_elem(3,2, 10.);
  m_src.copy_upper_triangle_into(m);
}


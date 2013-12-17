use sparse_matrix::*;
use lapack::lapack_int;
use lapack;

use std::vec;

#[test]
fn test_do_lapack_init() {
  lapack::init(); // TODO: Do this somewhere else, where it's gauranteed to be run before other tests as part of each test setup.
}


#[test]
fn test_3x4_push_and_get() {
  let mut m = SparseMatrix::new_with_capacities(12, 4, General);
  // row 0
  m.push(0,0, 0.);
  m.push(0,1, 1.);

  // row 1
  m.push(1,0, 3.);
  m.push(1,2, 5.);

  // row 2
  m.push(2,0, 6.);

  // row 3
  m.push(3,1, 10.);
  m.push(3,2, 11.);

  assert_eq!(m.get(0,0), 0.);
  assert_eq!(m.get(0,1), 1.);
  assert_eq!(m.get(0,2), 0.);
  
  assert_eq!(m.get(1,0), 3.);
  assert_eq!(m.get(1,1), 0.);
  assert_eq!(m.get(1,2), 5.);
  
  assert_eq!(m.get(2,0), 6.);
  assert_eq!(m.get(2,1), 0.);
  assert_eq!(m.get(2,2), 0.);
  
  assert_eq!(m.get(3,0), 0.);
  assert_eq!(m.get(3,1), 10.);
  assert_eq!(m.get(3,2), 11.);
}

#[test]
#[should_fail]
fn test_bad_values_capacity() {
  //      1 2 3
  // A =  2 2 0
  //      3 0 3 
  // b = [3 2 1]^t
  // sol = [0 1 1/3]^t
  let mut A = SparseMatrix::new_with_capacities(6, 3, General);
  A.push(0,0, 1.);
  A.push(0,1, 2.);
  A.push(0,2, 3.);
  A.push(1,0, 2.);
  A.push(1,1, 2.);
  A.push(2,0, 3.);
  A.push(2,2, 3.);
}

#[test]
fn test_3x4_csr3() {
  let mut m = SparseMatrix::new_with_capacities(12, 4, General);
  // row 0
  m.push(0,0, 0.);
  m.push(0,1, 1.);

  // row 1
  m.push(1,0, 3.);
  m.push(1,2, 5.);

  // row 2
  m.push(2,0, 6.);

  // row 3
  m.push(3,1, 10.);
  m.push(3,2, 11.);

  unsafe {
    let (vals, row_begins, cols) = m.csr3_ptrs();
    
    let vals_vec = vec::from_buf(vals, 7);
    assert_eq!(vals_vec, ~[0.,1.,3.,5.,6.,10.,11.]);

    let cols_vec = vec::from_buf(cols, 7);
    assert_eq!(cols_vec, ~[0,1,0,2,0,1,2 as lapack_int]);

    let row_begins_vec = vec::from_buf(row_begins, 5);
    assert_eq!(row_begins_vec, ~[0,2,4,5,7 as lapack_int]);
  } 

  assert_eq!(m.num_rows(), 4);
  assert_eq!(m.num_values(), 7);
}


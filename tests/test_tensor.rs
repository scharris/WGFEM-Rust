use tensor::*;
use common::R;

#[test]
fn test_tensor3_constr_from_elem() {
  let t = Tensor3::from_elem(3,4,5, 99.);
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        assert_eq!(t.get(i,j,k), 99.);
      }
    }
  }
}

#[test]
fn test_tensor3_constr_from_fn() {
  let t = Tensor3::from_fn(3,4,5, |i,j,k| 1000.*(i as R) + 100.*(j as R) + (k as R));
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        assert_eq!(t.get(i,j,k), 1000.*(i as R) + 100.*(j as R) + (k as R));
      }
    }
  }
}

#[test]
fn test_tensor4_constr_from_elem() {
  let t = Tensor4::from_elem(3,4,5,6, 99.);
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        for l in range(0u,6) {
          assert_eq!(t.get(i,j,k,l), 99.);
        }
      }
    }
  }
}

#[test]
fn test_tensor4_constr_from_fn() {
  let t = Tensor4::from_fn(3,4,5,6, |i,j,k,l| 1000.*(i as R) + 100.*(j as R) + 10.*(k as R) + (l as R));
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        for l in range(0u,6) {
          assert_eq!(t.get(i,j,k,l), 1000.*(i as R) + 100.*(j as R) + 10.*(k as R) + (l as R));
        }
      }
    }
  }
}

#[test]
fn test_tensor5_constr_from_elem() {
  let t = Tensor5::from_elem(3,4,5,6,7, 99.);
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        for l in range(0u,6) {
          for m in range(0u,7) {
            assert_eq!(t.get(i,j,k,l,m), 99.);
          }
        }
      }
    }
  }
}

#[test]
fn test_tensor5_constr_from_fn() {
  let t = Tensor5::from_fn(3,4,5,6,7, |i,j,k,l,m| 1000.*(i as R) + 100.*(j as R) + 10.*(k as R) + (l as R) + (m as R)/10.);
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        for l in range(0u,6) {
          for m in range(0u,7) {
            assert_eq!(t.get(i,j,k,l,m), 1000.*(i as R) + 100.*(j as R) + 10.*(k as R) + (l as R) + (m as R)/10.);
          }
        }
      }
    }
  }
}


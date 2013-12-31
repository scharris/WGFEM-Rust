use storage_by_ints::*;
use common::R;

#[test]
fn test_StorageByInts2_constr_from_fn() {
  let t = StorageByInts2::from_fn(3,4, |i,j| 1000.*(i as R) + 100.*(j as R));
  for i in range(0u,3) {
    for j in range(0u,4) {
      assert_eq!(t.get(i,j), 1000.*(i as R) + 100.*(j as R));
    }
  }
}

#[test]
fn test_StorageByInts2_set() {
  let mut t = StorageByInts2::from_fn(3,4, |i,j| 1000.*(i as R) + 100.*(j as R));
  t.set(1,1, -99.);
  for i in range(0u,3) {
    for j in range(0u,4) {
      if i == 1 && j == 1 {
        assert_eq!(t.get(i,j), -99.);
      } else {
        assert_eq!(t.get(i,j), 1000.*(i as R) + 100.*(j as R));
      }
    }
  }
}

#[test]
fn test_StorageByInts2_constr_from_elem() {
  let mut t = StorageByInts2::from_elem(3,4, 99.);
  for i in range(0u,3) {
    for j in range(0u,4) {
      assert_eq!(t.get(i,j), 99.);
    }
  }
}


#[test]
fn test_StorageByInts3_constr_from_fn() {
  let t = StorageByInts3::from_fn(3,4,5, |i,j,k| 1000.*(i as R) + 100.*(j as R) + (k as R));
  for i in range(0u,3) {
    for j in range(0u,4) {
      for k in range(0u,5) {
        assert_eq!(t.get(i,j,k), 1000.*(i as R) + 100.*(j as R) + (k as R));
      }
    }
  }
}

#[test]
fn test_StorageByInts4_constr_from_fn() {
  let t = StorageByInts4::from_fn(3,4,5,6, |i,j,k,l| 1000.*(i as R) + 100.*(j as R) + 10.*(k as R) + (l as R));
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
fn test_StorageByInts5_constr_from_fn() {
  let t = StorageByInts5::from_fn(3,4,5,6,7, |i,j,k,l,m| 1000.*(i as R) + 100.*(j as R) + 10.*(k as R) + (l as R) + (m as R)/10.);
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


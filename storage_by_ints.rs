use std::vec;

struct StorageByInts2<T> {
  data: ~[T],
  size_0: uint,
  size_1: uint,
}

impl <T:Clone> StorageByInts2<T> {

  pub fn from_elem(size_0: uint, size_1: uint, elem: T) -> StorageByInts2<T> {
    let sz = size_0 * size_1;
    let mut data = vec::with_capacity(sz);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        data.push(elem.clone());
      }
    }
    StorageByInts2 {
      data: data,
      size_0: size_0,
      size_1: size_1,
    }
  }

  pub fn from_fn(size_0: uint, size_1: uint,
                 f: |i0:uint, i1:uint| -> T) -> StorageByInts2<T> {
    let mut data = vec::with_capacity(size_0 * size_1);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        data.push(f(i0, i1));
      }
    }
    StorageByInts2 {
      data: data,
      size_0: size_0,
      size_1: size_1,
    }
  }

  #[inline]
  pub fn get(&self, i0: uint, i1: uint) -> T {
    self.data[i0 * self.size_1 + i1].clone()
  }
  
  #[inline]
  pub fn set(&mut self, i0: uint, i1: uint, val: T) {
    self.data[i0 * self.size_1 + i1] = val;
  }
  
}

struct StorageByInts3<T> {
  data: ~[T],
  size_0: uint,
  size_1: uint,
  size_2: uint,
}

impl <T:Clone> StorageByInts3<T> {

  pub fn from_fn(size_0: uint, size_1: uint, size_2: uint,
                 f: |i0:uint, i1:uint, i2:uint| -> T) -> StorageByInts3<T> {
    let mut data = vec::with_capacity(size_0 * size_1 * size_2);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        for i2 in range(0, size_2) {
          data.push(f(i0, i1, i2));
        }
      }
    }
    StorageByInts3 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
    }
  }

  #[inline]
  pub fn get(&self, i0: uint, i1: uint, i2: uint) -> T {
    self.data[(i0 * self.size_1 + i1)*self.size_2 + i2].clone()
  }
}

struct StorageByInts4<T> {
  data: ~[T],
  size_0: uint,
  size_1: uint,
  size_2: uint,
  size_3: uint,
}

impl <T:Clone> StorageByInts4<T> {

  pub fn from_elem(size_0: uint, size_1: uint, size_2: uint, size_3: uint, elem: T) -> StorageByInts4<T> {
    let sz = size_0 * size_1 * size_2 * size_3;
    let mut data = vec::with_capacity(sz);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        for i2 in range(0, size_2) {
          for i3 in range(0, size_3) {
            data.push(elem.clone());
          }
        }
      }
    }
    StorageByInts4 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
    }
  }

  pub fn from_fn(size_0: uint, size_1: uint, size_2: uint, size_3: uint,
                 f: |i0:uint, i1:uint, i2:uint, i3: uint| -> T) -> StorageByInts4<T> {
    let sz = size_0 * size_1 * size_2 * size_3;
    let mut data = vec::with_capacity(sz);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        for i2 in range(0, size_2) {
          for i3 in range(0, size_3) {
            data.push(f(i0, i1, i2, i3));
          }
        }
      }
    }
    StorageByInts4 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
    }
  }

  #[inline]
  pub fn get(&self, i0: uint, i1: uint, i2: uint, i3: uint) -> T {
    self.data[((i0 * self.size_1 + i1) * self.size_2 + i2) * self.size_3 + i3].clone()
  }
}


struct StorageByInts5<T> {
  data: ~[T],
  size_0: uint,
  size_1: uint,
  size_2: uint,
  size_3: uint,
  size_4: uint,
}

impl <T:Clone> StorageByInts5<T> {

  pub fn from_fn(size_0: uint, size_1: uint, size_2: uint, size_3: uint, size_4: uint,
                 f: |i0:uint, i1:uint, i2:uint, i3: uint, i4: uint| -> T) -> StorageByInts5<T> {
    let sz = size_0 * size_1 * size_2 * size_3 * size_4;
    let mut data = vec::with_capacity(sz);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        for i2 in range(0, size_2) {
          for i3 in range(0, size_3) {
            for i4 in range(0, size_4) {
              data.push(f(i0, i1, i2, i3, i4));
            }
          }
        }
      }
    }
    StorageByInts5 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
      size_4: size_4,
    }
  }

  #[inline]
  pub fn get(&self, i0: uint, i1: uint, i2: uint, i3: uint, i4: uint) -> T {
    self.data[(((i0 * self.size_1 + i1) * self.size_2 + i2) * self.size_3 + i3) * self.size_4 + i4].clone()
  }
}


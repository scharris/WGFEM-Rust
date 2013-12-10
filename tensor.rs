use common::{R};

use std::vec;


struct Tensor3 {
  data: ~[R],
  size_0: uint,
  size_1: uint,
  size_2: uint,
}

impl Tensor3 {

/*
  pub fn from_elem(size_0: uint, size_1: uint, size_2: uint, elem: R) -> Tensor3 {
    let sz = size_0 * size_1 * size_2;
    Tensor3 {
      data: vec::from_elem(sz, elem),
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
    }
  }
*/ 

  pub fn from_fn(size_0: uint, size_1: uint, size_2: uint, f: |i0:uint, i1:uint, i2:uint| -> R) -> Tensor3 {
    let mut data = vec::with_capacity(size_0 * size_1 * size_2);
    for i0 in range(0, size_0) {
      for i1 in range(0, size_1) {
        for i2 in range(0, size_2) {
          data.push(f(i0, i1, i2));
        }
      }
    }
    Tensor3 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
    }
  }

/*
  #[inline]
  pub fn set(&mut self, i0: uint, i1: uint, i2: uint, val: R) {
    self.data[(i0 * self.size_1 + i1)*self.size_2 + i2] = val;
  }
*/
  
  #[inline]
  pub fn get(&self, i0: uint, i1: uint, i2: uint) -> R {
    self.data[(i0 * self.size_1 + i1)*self.size_2 + i2]
  }
}

struct Tensor4 {
  data: ~[R],
  size_0: uint,
  size_1: uint,
  size_2: uint,
  size_3: uint,
}

impl Tensor4 {

  pub fn from_elem(size_0: uint, size_1: uint, size_2: uint, size_3: uint, elem: R) -> Tensor4 {
    let sz = size_0 * size_1 * size_2 * size_3;
    Tensor4 {
      data: vec::from_elem(sz, elem),
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
    }
  }
  
  pub fn from_fn(size_0: uint, size_1: uint, size_2: uint, size_3: uint, f: |i0:uint, i1:uint, i2:uint, i3: uint| -> R) -> Tensor4 {
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
    Tensor4 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
    }
  }

/*
  #[inline]
  pub fn set(&mut self, i0: uint, i1: uint, i2: uint, i3: uint, val: R) {
    self.data[((i0 * self.size_1 + i1) * self.size_2 + i2) * self.size_3 + i3] = val;
  }
*/  

  #[inline]
  pub fn get(&self, i0: uint, i1: uint, i2: uint, i3: uint) -> R {
    self.data[((i0 * self.size_1 + i1) * self.size_2 + i2) * self.size_3 + i3]
  }
}


struct Tensor5 {
  data: ~[R],
  size_0: uint,
  size_1: uint,
  size_2: uint,
  size_3: uint,
  size_4: uint,
}

impl Tensor5 {

/*
  pub fn from_elem(size_0: uint, size_1: uint, size_2: uint, size_3: uint, size_4: uint, elem: R) -> Tensor5 {
    let sz = size_0 * size_1 * size_2 * size_3 * size_4;
    Tensor5 {
      data: vec::from_elem(sz, elem),
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
      size_4: size_4,
    }
  }
*/

  pub fn from_fn(size_0: uint, size_1: uint, size_2: uint, size_3: uint, size_4: uint, f: |i0:uint, i1:uint, i2:uint, i3: uint, i4: uint| -> R) -> Tensor5 {
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
    Tensor5 {
      data: data,
      size_0: size_0,
      size_1: size_1,
      size_2: size_2,
      size_3: size_3,
      size_4: size_4,
    }
  }

/*
  #[inline]
  pub fn set(&mut self, i0: uint, i1: uint, i2: uint, i3: uint, i4: uint, val: R) {
    self.data[(((i0 * self.size_1 + i1) * self.size_2 + i2) * self.size_3 + i3) * self.size_4 + i4] = val;
  }
*/
  
  #[inline]
  pub fn get(&self, i0: uint, i1: uint, i2: uint, i3: uint, i4: uint) -> R {
    self.data[(((i0 * self.size_1 + i1) * self.size_2 + i2) * self.size_3 + i3) * self.size_4 + i4]
  }
}


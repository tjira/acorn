use std::fmt;
use std::ops;

use number::Number;

#[derive(Debug, Default)]
pub struct Tensor<T> {
    data: Vec<T>, shape: Vec<usize>,
}

impl<T: Number> Tensor<T> {
    pub fn init(shape: Vec<usize>) -> Self {
        Tensor{data: vec![T::default(); shape.iter().product()], shape}
    }
    pub fn linspace(start: T, end: T, shape: Vec<usize>) -> Self {
        let mut data = Vec::with_capacity(shape.iter().product());
        for i in 0..data.capacity() {
            data.push(start + (end - start) * (i as f64 / (data.capacity() - 1) as f64));
        }
        Tensor{data, shape: shape}
    }
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl<T: Number> ops::Add for Tensor<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Tensor{data: self.data.iter().zip(rhs.data.iter()).map(|(&x, &y)| x + y).collect(), shape: self.shape}
    }
}

impl<T: fmt::Display> fmt::Display for Tensor<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[")?;
        for i in 0..self.data.len() {
            write!(f, "{}{}", if i > 0 {" "} else {""}, self.data[i])?;
        }
        write!(f, "]")
    }
}

impl<T> ops::Index<usize> for Tensor<T> {
    type Output = T;
    fn index(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<T> ops::IndexMut<usize> for Tensor<T> {
    fn index_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }
}

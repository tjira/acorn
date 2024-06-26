use std::ops; use number::Real;

#[derive(Clone, Debug, Default)]
pub struct Array<T> {
    data: Vec<T>, shape: Vec<usize>
}

impl<T: ops::Add<Output=T> + ops::Mul<Output=T> + Clone + Copy + Default> Array<T> {
    // functions that create new data
    pub fn init(shape: Vec<usize>) -> Array<T> {
        Array{data: vec![T::default(); shape.iter().product()], shape: shape}
    }
    pub fn new(data: Vec<T>) -> Array<T> {
        Array{shape: vec![data.len()], data: data}
    }

    // immutable references to the data
    pub fn shape(&self) -> &Vec<usize> {
        &self.shape
    }
    pub fn size(&self) -> usize {
        self.data.len()
    }

    // mutable references to the data
    pub fn data(&mut self) -> &mut Vec<T> {
        &mut self.data
    }

    // vector functions
    pub fn apply<F: Fn(T) -> T>(self, f: F) -> Array<T> {
        Array{data: self.data.into_iter().map(f).collect(), shape: self.shape}
    }
    pub fn cast<U: From<T>>(self) -> Array<U> {
        Array{data: self.data.into_iter().map(|x| x.into()).collect(), shape: self.shape}
    }
    pub fn conj(self) -> Array<T> where T: num::complex::ComplexFloat {
        Array{data: self.data.into_iter().map(|x| x.conj()).collect(), shape: self.shape}
    }
    pub fn dot(&self, rhs: &Array<T>) -> T {
        self.data.iter().copied().zip(rhs.data.iter().copied()).map(|(x, y)| x * y).fold(T::default(), |a, x| a + x)
    }
    pub fn sum(&self) -> T where T: From<Real> {
        self.data.iter().copied().fold(0.0.into(), ops::Add::add)
    }
}

impl<T: ops::Mul<Output=T>> ops::Mul<Array<T>> for Array<T> {
    type Output = Array<T>;
    fn mul(self, rhs: Array<T>) -> Array<T> {
        Array{data: self.data.into_iter().zip(rhs.data.into_iter()).map(|(x, y)| x * y).collect(), shape: self.shape}
    }
}

impl<T: ops::Mul<Output=T> + From<Real> + Copy> ops::Mul<Array<T>> for Real {
    type Output = Array<T>;
    fn mul(self, rhs: Array<T>) -> Array<T> {
        Array{data: rhs.data.into_iter().map(|x| x * self.into()).collect(), shape: rhs.shape}
    }
}

impl<T: ops::Mul<Output=T> + From<Real>> ops::Mul<Real> for Array<T> {
    type Output = Array<T>;
    fn mul(self, rhs: Real) -> Array<T> {
        Array{data: self.data.into_iter().map(|x| x * rhs.into()).collect(), shape: self.shape}
    }
}

impl<T> ops::Index<usize> for Array<T> {
    type Output = T;
    fn index(&self, index: usize) -> &T {
        &self.data[index]
    }
}

use std::fmt; use std::ops; use number::{Real, RealTrait, ComplexTrait, NumberTrait};

#[derive(Clone, Debug, Default)]
pub struct Array<T> {
    data: Vec<T>, shape: Vec<usize>
}

impl<T: NumberTrait> Array<T> {
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
    pub fn cast<U: NumberTrait + From<T>>(self) -> Array<U> {
        Array{data: self.data.into_iter().map(|x| x.into()).collect(), shape: self.shape}
    }
    pub fn conj(self) -> Array<T> where T: ComplexTrait {
        Array{data: self.data.into_iter().map(|x| x.conj()).collect(), shape: self.shape}
    }
    pub fn meshgrided(self, dim: u32) -> Array<T> {
        Array{data: (0..self.size().pow(dim)).flat_map(|_| self.data.iter().cloned()).collect(), shape: (0..dim).flat_map(|_| self.shape.iter().cloned()).collect()}

    }
    pub fn reshaped(self, shape: Vec<usize>) -> Array<T> {
        Array{data: self.data, shape: shape}
    }
    pub fn sqrt(self) -> Array<T> where T: RealTrait {
        Array{data: self.data.into_iter().map(|x| x.sqrt()).collect(), shape: self.shape}
    }
    pub fn dot(&self, rhs: &Array<T>) -> T {
        self.data.iter().copied().zip(rhs.data.iter().copied()).map(|(x, y)| x * y).fold(T::default(), |a, x| a + x)
    }
    pub fn sum(&self) -> T {
        self.data.iter().copied().fold(T::default(), ops::Add::add)
    }
}

impl<T: NumberTrait> ops::Mul<Array<T>> for Array<T> {
    type Output = Array<T>;
    fn mul(self, rhs: Array<T>) -> Array<T> {
        Array{data: self.data.into_iter().zip(rhs.data.into_iter()).map(|(x, y)| x * y).collect(), shape: self.shape}
    }
}

impl<T: NumberTrait + From<Real>> ops::Mul<Array<T>> for Real {
    type Output = Array<T>;
    fn mul(self, rhs: Array<T>) -> Array<T> {
        Array{data: rhs.data.iter().copied().map(|x| <Real as Into<T>>::into(self) * x).collect(), shape: rhs.shape}
    }
}

impl<T: NumberTrait + From<Real>> ops::Mul<Real> for Array<T> {
    type Output = Array<T>;
    fn mul(self, rhs: Real) -> Array<T> {
        Array{data: self.data.iter().copied().map(|x| x * <Real as Into<T>>::into(rhs)).collect(), shape: self.shape}
    }
}

impl<T: NumberTrait> ops::Index<usize> for Array<T> {
    type Output = T;
    fn index(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<T: NumberTrait + fmt::Display> fmt::Display for Array<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..self.size() {
            write!(f, "{:>20.14}\n", self.data[i])?;
        } Ok(())
    }
}

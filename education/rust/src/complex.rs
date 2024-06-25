use std::fmt;
use std::ops;

use number::Number;

#[derive(Clone, Copy, Debug, Default)]
pub struct Complex<T> {
    real: T, imag: T,
}

impl<T: Number + From<f64> + Into<f64>> Complex<T> {
    pub fn abs(&self) -> T {
        <f64 as Into<T>>::into((self.real * self.real + self.imag * self.imag).into().sqrt())
    }
    pub fn conj(&self) -> Self {
        Complex{real: self.real, imag: -self.imag}
    }
    pub fn exp(&self) -> Self {
        Complex{real: (self.real.into().exp() * self.imag.into().cos()).into(), imag: (self.real.into().exp() * self.imag.into().sin()).into()}
    }
    pub fn i() -> Self {
        Complex{real: T::default(), imag: T::default() + <f64 as Into<T>>::into(1_f64)}
    }
    pub fn new(real: T, imag: T) -> Self {
        Complex{real, imag}
    }
}

impl<T: Number> ops::Add for Complex<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Complex{real: self.real + rhs.real, imag: self.imag + rhs.imag}
    }
}

impl<T: Number + From<f64>> ops::Add<f64> for Complex<T> {
    type Output = Self;
    fn add(self, rhs: f64) -> Self {
        Complex{real: self.real + <f64 as Into<T>>::into(rhs), imag: self.imag}
    }
}

impl<T: Number + From<f64>> ops::Add<Complex<T>> for f64 {
    type Output = Complex<T>;
    fn add(self, rhs: Complex<T>) -> Complex<T> {
        Complex{real: <f64 as Into<T>>::into(self) + rhs.real, imag: rhs.imag}
    }
}

impl<T: fmt::Display> fmt::Display for Complex<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}+{}i", self.real, self.imag)
    }
}


impl<T: Number> ops::Div<f64> for Complex<T> {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Complex{real: self.real / rhs, imag: self.imag / rhs}
    }
}

impl<T: Number> ops::Mul for Complex<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Complex{real: self.real * rhs.real - self.imag * rhs.imag, imag: self.real * rhs.imag + self.imag * rhs.real}
    }
}

impl<T: Number> ops::Mul<f64> for Complex<T> {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Complex{real: self.real * rhs, imag: self.imag * rhs}
    }
}

impl<T: Number + From<f64>> ops::Mul<Complex<T>> for f64 {
    type Output = Complex<T>;
    fn mul(self, rhs: Complex<T>) -> Complex<T> {
        Complex{real: <f64 as Into<T>>::into(self) * rhs.real, imag: <f64 as Into<T>>::into(self) * rhs.imag}
    }
}

impl<T: Number> ops::Sub for Complex<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Complex{real: self.real - rhs.real, imag: self.imag - rhs.imag}
    }
}

impl<T: Number + From<f64>> ops::Sub<f64> for Complex<T> {
    type Output = Self;
    fn sub(self, rhs: f64) -> Self {
        Complex{real: self.real - <f64 as Into<T>>::into(rhs), imag: self.imag}
    }
}

impl<T: Number + From<f64>> ops::Sub<Complex<T>> for f64 {
    type Output = Complex<T>;
    fn sub(self, rhs: Complex<T>) -> Complex<T> {
        Complex{real: <f64 as Into<T>>::into(self) - rhs.real, imag: -rhs.imag}
    }
}

impl<T: Number> ops::Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Complex{real: -self.real, imag: -self.imag}
    }
}

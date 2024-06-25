use std::ops;

use complex::Complex;

pub trait Number: Clone + Copy + Default + ops::Add<Output = Self> + ops::Mul<Output = Self> + ops::Sub<Output = Self> + ops::Add<f64, Output = Self> + ops::Sub<f64, Output = Self> + ops::Mul<f64, Output = Self> + ops::Div<f64, Output = Self> + ops::Neg<Output = Self> {}

impl Number for f64 {} impl<T: Number + From<f64>> Number for Complex<T> {}

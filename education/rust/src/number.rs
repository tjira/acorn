// define the used types
pub type Real = f64; pub type Complex = num::complex::Complex<Real>;

// trait for all allowed real and complex numbers
pub trait NumberTrait: Copy + Default + num::traits::Num {}
pub trait RealTrait: NumberTrait + num::traits::Float {}
pub trait ComplexTrait: NumberTrait + num::complex::ComplexFloat {}

// implement traits for reals
impl NumberTrait for Real {}
impl RealTrait for Real {}

// implement traits for complex
impl NumberTrait for Complex {}
impl ComplexTrait for Complex {}

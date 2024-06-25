extern crate fftw; use fourier::fftw::plan::C2CPlan;

use number::Complex; use vector::Vector;

pub fn fftcc(mut data: Vector<Complex>) -> Vector<Complex> {
    // define the plan of the fft
    let mut plan: fftw::plan::C2CPlan64 = C2CPlan::aligned(&data.shape(), fftw::types::Sign::Forward, fftw::types::Flag::MEASURE).unwrap();
    
    // create a new vector to store the output
    let mut output = Vector::init(data.shape().clone());

    // perform the fft and return the output in the correct format
    plan.c2c(&mut data.data(), &mut output.data()).unwrap(); return output;
}

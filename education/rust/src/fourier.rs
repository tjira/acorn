extern crate fftw; use fourier::fftw::plan::C2CPlan;

use number::Real; use number::Complex; use array::Array;

pub fn fftcc(mut data: Array<Complex>) -> Array<Complex> {
    // define the plan of the fft
    let mut plan: fftw::plan::C2CPlan64 = C2CPlan::aligned(&data.shape(), fftw::types::Sign::Forward, fftw::types::Flag::MEASURE).unwrap();
    
    // create a new vector to store the output
    let mut output = Array::init(data.shape().clone());

    // perform the fft and return the output in the correct format
    plan.c2c(&mut data.data(), &mut output.data()).unwrap(); return output;
}

pub fn ifftcc(mut data: Array<Complex>) -> Array<Complex> {
    // define the plan of the ifft
    let mut plan: fftw::plan::C2CPlan64 = C2CPlan::aligned(&data.shape(), fftw::types::Sign::Backward, fftw::types::Flag::MEASURE).unwrap();
    
    // create a new vector to store the output
    let mut output = Array::init(data.shape().clone());

    // perform the ifft and return the output in the correct format
    plan.c2c(&mut data.data(), &mut output.data()).unwrap(); return 1.0 / output.size() as Real * output;
}

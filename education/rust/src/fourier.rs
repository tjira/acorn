// extern crate imports
extern crate fftw;

// crate includes
use num::complex::Complex; use ndarray::{Array, Ix1}; use fourier::fftw::plan::C2CPlan;

pub fn fftcc(mut data: Array<Complex<f64>, Ix1>, d: usize) -> Array<Complex<f64>, Ix1> {
    // calculate the shape
    let shape = vec![(data.len() as f64).powf(1.0 / d as f64).round() as usize; d];

    // define the plan of the fft
    let mut plan: fftw::plan::C2CPlan64 = C2CPlan::aligned(&shape, fftw::types::Sign::Forward, fftw::types::Flag::MEASURE).unwrap();
    
    // create a new vector to store the output
    let mut output = Array::<Complex<f64>, Ix1>::zeros(data.len());

    // perform the fft and return the output in the correct format
    plan.c2c(data.as_slice_mut().unwrap(), output.as_slice_mut().unwrap()).unwrap(); output
}

pub fn ifftcc(mut data: Array<Complex<f64>, Ix1>, d: usize) -> Array<Complex<f64>, Ix1> {
    // calculate the shape
    let shape = vec![(data.len() as f64).powf(1.0 / d as f64).round() as usize; d];

    // define the plan of the ifft
    let mut plan: fftw::plan::C2CPlan64 = C2CPlan::aligned(&shape, fftw::types::Sign::Backward, fftw::types::Flag::MEASURE).unwrap();
    
    // create a new vector to store the output
    let mut output = Array::<Complex<f64>, Ix1>::zeros(data.len());

    // perform the ifft and return the output in the correct format
    plan.c2c(data.as_slice_mut().unwrap(), output.as_slice_mut().unwrap()).unwrap(); output / data.len() as f64
}

pub fn fftfreq(size: usize, step: f64) -> Array<f64, Ix1> {
    // create the container for the data
    let mut data = Array::<f64, Ix1>::from_elem(size, 2.0 * std::f64::consts::PI / (step * size as f64));

    // fill the container with the data
    for i in 0..size {
        data[i] = data[i] * (i as f64 - if i < size / 2 {0.0} else {size as f64})
    }

    data // return the data
}

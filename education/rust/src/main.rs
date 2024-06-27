// extern crate imports
extern crate num; extern crate ndarray;

// module imports
mod fourier; mod wavefunction;

// system imports
use std::time::{Instant};

// crate function and object imports
use num::complex::Complex; use ndarray::{Array, ArrayView, Axis, Ix1, Ix2};

// object imports
use wavefunction::{State, Wavefunction};

// function imports
use fourier::{fftfreq};

fn main() {
    // start the timer
    let start = Instant::now();

    // define the variables
    let d: usize = 2; let iters: usize = 1000; let points: usize = 64; let rmin: f64 = -8.0; let rmax: f64 = 8.0;

    // define potential and initial wavefunction functions
    let uf = |r: ArrayView<f64, Ix1>| 0.5 * r.mapv(|x| x * x).sum(); let wf = |r: ArrayView<f64, Ix1>| (-r.mapv(|x| x * x).sum()).exp();

    // define the r-space and k-space array in one dimension
    let rr = Array::<f64, Ix1>::linspace(rmin, rmax, points); let kk = fftfreq(rr.len(), rr[1] - rr[0]);

    // define the r-space and k-space array in d dimensions
    let mut r = Array::<f64, Ix2>::zeros((points.pow(d as u32), d)); let mut k = Array::<f64, Ix2>::zeros((points.pow(d as u32), d));

    // fill the r-space and k-space arrays
    for i in 0..r.shape()[0] {
        for j in 0..r.shape()[1] {
            r[[i, d - j - 1]] = rr[i / points.pow(j as u32) % points];
            k[[i, d - j - 1]] = kk[i / points.pow(j as u32) % points];
        }
    }

    // define the k-space array squared
    let ksq = k.axis_iter(Axis(0)).map(|x| Complex::from(x.mapv(|x| x * x).sum())).collect::<Array::<Complex<f64>, Ix1>>();

    // define the potential array
    let u = r.axis_iter(Axis(0)).map(uf).collect::<Array<f64, Ix1>>().mapv(|x| Complex::<f64>::from(x));

    // define the initial wavefunction
    let mut w = State::new(&r, wf, 1.0).normalized();

    // define the propagators
    let (rp, kp) = w.propagators(&ksq, &u, 0.1);

    // propagation loop
    for _ in 0..iters {
        w = w.propagated(&rp, &kp).normalized();
    }

    // print the energy
    println!("FINAL WFN ENERGY: {:?}", Wavefunction::new(vec![w]).hamilton(&ksq, &u).re);

    // print the elapsed time
    println!("ELAPSED TIME: {:?}", start.elapsed());
}

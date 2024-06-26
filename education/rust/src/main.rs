mod array; mod fourier; mod number; mod vector; mod wavefunction;

// system
use std::time::{Instant};

// types
use number::Complex; use array::Array; use wavefunction::Wavefunction;

// functions
use vector::{fftfreq, linspace};

fn main() {
    // start the timer
    let start = Instant::now();

    // define the variables
    let dim: u32 = 1; let iters = 1000; let points = 4; let rmin = -8.0; let rmax = 8.0;

    // define the r-space and k-space arrays
    let r = Array::new(linspace(rmin, rmax, points)); let ksq = Array::new(fftfreq(r.size(), r[1] - r[0])).apply(|x| x * x).cast::<Complex>();

    // define the potential array
    let u = r.clone().apply(|x| 0.5 * x * x).cast::<Complex>();

    // define the initial wavefunction
    let mut w = Wavefunction::new(&r, |x| (-x * x).exp(), 1.0).normalized();

    // define the propagators
    let (rp, kp) = w.propagators(&ksq, &u, 0.1);

    // propagation loop
    for _ in 0..iters {

        // propagate the wavefunction
        w = w.propagated(&rp, &kp).normalized();

    }

    // print the energy
    println!("{:?}", w.energy(&ksq, &u).norm());

    // print the elapsed time
    println!("{:?}", start.elapsed())
}

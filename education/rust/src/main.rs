// extern crate imports
extern crate num; extern crate ndarray;

// module imports
mod fourier; mod matrix; mod qdyn; mod wavefunction;

// system includes
use std::time::Instant;

// crate includes
use ndarray::{ArrayView, Ix1};

// personal includes
use matrix::writemat;

fn main() {
    // define the qdyn option struct
    let qdyn = qdyn::QuantumDynamics{d: 2, iters: 1000, nstate: 9, points: 64, rmin: -8.0, rmax: 8.0, imaginary: true, savewfn: true};

    // define potential and initial wavefunction functions
    let ufv = vec![|r: ArrayView<f64, Ix1>| 0.5 * r.mapv(|x| x * x).sum()]; let wfv = vec![|r: ArrayView<f64, Ix1>| (-r.mapv(|x| (x - 1.0) * (x - 1.0)).sum()).exp()];

    // start the timer
    let start = Instant::now();

    // perform the dynamics and save the potential
    let (u, wt, ei) = qdyn.run(ufv, wfv); writemat("U.mat", &u);

    // print the energy and save the wavefunctions
    for (i, e) in ei.iter().enumerate() {
        println!("FINAL WFN {:02} ENERGY: {:.14}", i, e);
    }

    // save the wavefunctions
    if wt.is_some() {
        for (i, w) in wt.unwrap().iter().enumerate() {
            writemat(&format!("PSI_DIA_{:02}.mat", i), &w)
        }
    }

    // print the elapsed time
    println!("\nELAPSED TIME: {:?}", start.elapsed())
}

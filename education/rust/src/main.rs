extern crate ndarray; extern crate num_complex; use num_complex::{Complex64 as c64}; use ndarray::prelude::{array, s}; use ndarray::{Array1 as Vector, ArrayView1 as VectorView, Array2 as Matrix, Axis, concatenate};

mod qdyn;

use std::fs; use std::io::{BufWriter, Write}; use std::ops; use std::time;

pub fn fftfreq(size: usize, step: f64) -> Vector<f64> {
    (0..size).map(|i| 2.0 * std::f64::consts::PI / (step * size as f64) * (i as f64 - if i < size / 2 {0.0} else {size as f64})).collect()
}

pub fn writemat(path: &str, mat: &Matrix<f64>) {
    // open the file
    let file = fs::File::create(path).expect("UNABLE TO WRITE");

    // create the buffer
    let mut buffer = BufWriter::new(file);

    // write the header
    write!(buffer, "{} {}\n", mat.nrows(), mat.ncols()).expect("UNABLE TO WRITE");

    for i in 0..mat.nrows() {
        for j in 0..mat.ncols() {
            write!(buffer, "{:20.14}{}", mat[[i, j]], if j == mat.ncols() - 1 {"\n"} else {" "}).expect("UNABLE TO WRITE");
        }
    }
}

fn main() {
    // define the quantum dynamics parameters
    let qdyn = qdyn::QuantumDynamics{d: 1, iters: 1000, nstate: 9, points: 8, dt: 0.1, mass: 1.0, rmax: 8.0, rmin: -8.0, imaginary: true, savewfn: true};

    // define the potential function
    let uf = |r: &VectorView<f64>| -> Matrix<f64> {array![
        [0.5 * r.mapv(|x| x.powi(2)).sum()]
    ]};

    // define the wavefunction function
    let wf = |r: &Vector<f64>| -> Vector<c64> {array!
        [(-0.5 * r.mapv(|x| x.powi(2)).sum()).exp()].mapv(|x| c64::new(x, 0.0))
    };

    // println!("{:?}", wf(&array![2.0]));

    // define the r-space and k-space grids
    let mut r = Matrix::<f64>::zeros((qdyn.points.pow(qdyn.d as u32), qdyn.d));
    let mut k = Matrix::<f64>::zeros((qdyn.points.pow(qdyn.d as u32), qdyn.d));

    // define the delta r
    let dr = (qdyn.rmax - qdyn.rmin) / (qdyn.points - 1) as f64;

    // assign the linspace values to the 1D sector of the r-space grid
    r.slice_mut(s![0..qdyn.points, qdyn.d - 1]).assign(&Vector::linspace(qdyn.rmin, qdyn.rmax, qdyn.points));

    // assign the linspace values to the 1D sector of the k-space grid
    k.slice_mut(s![0..qdyn.points, qdyn.d - 1]).assign(&fftfreq(qdyn.points, dr));

    // fill the additional dimensions
    for i in 0..r.shape()[0] {
        for j in 0..r.shape()[1] {
            r[[i, qdyn.d - j - 1]] = r[[i / qdyn.points.pow(j as u32) % qdyn.points, qdyn.d - 1]];
            k[[i, qdyn.d - j - 1]] = k[[i / qdyn.points.pow(j as u32) % qdyn.points, qdyn.d - 1]];
        }
    }

    // define the wavefunction
    // let mut psi = Matrix::<c64>::zeros((qdyn.points.pow(qdyn.d as u32), 1));
    let mut w = r.axis_iter(Axis(0)).map(|x| wf(&x.to_owned())).collect::<Vec<Vector<c64>>>();

    // assign the wavefunction values
    // for i in 0..psi.shape()[0] {
    //     psi.slice_mut(s![i, ..]).assign(&wf(&r.slice(s![i, ..]).to_owned()));
    // }

    // writemat("PSI.mat", &concatenate![Axis(1), r, w.mapv(|x| x.re), w.mapv(|x| x.im)]);
    println!("{:?}", w);

}

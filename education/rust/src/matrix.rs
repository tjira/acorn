// extern crate imports
extern crate nalgebra; extern crate nshare;

// system includes
use std::io::{BufWriter, Write}; use std::fs::File;

// crate includes
use ndarray::{Array, Ix1, Ix2}; use matrix::nalgebra::SymmetricEigen; use matrix::nshare::{MutNdarray1, ToNalgebra, ToNdarray2};

pub fn eigh(m: &Array<f64, Ix2>) -> (Array<f64, Ix2>, Array<f64, Ix1>) {
    let eigen = SymmetricEigen::new(m.clone().into_nalgebra()); (eigen.eigenvectors.into_ndarray2(), eigen.eigenvalues.clone().mut_ndarray1().to_owned())
}

pub fn writemat(path: &str, mat: &Array<f64, Ix2>) {
    // open the file
    let file = File::create(path).expect("UNABLE TO WRITE FILE");

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

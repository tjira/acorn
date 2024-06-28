// system includes
use std::io::{BufWriter, Write}; use std::fs::File;

// crate includes
use ndarray::{Array, Ix2};

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

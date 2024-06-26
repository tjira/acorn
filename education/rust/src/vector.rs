use number::{Real, NumberTrait};

pub fn fftfreq(size: usize, step: Real) -> Vec<Real> {
    // create the container for the data
    let mut data = vec![2.0 * std::f64::consts::PI / (step * size as Real); size];

    // fill the container with the data
    for i in 0..size {
        data[i] = data[i] * (i as Real - if i < size / 2 {0.0} else {size as Real});
    }

    return data; // return the data
}

pub fn linspace(start: Real, end: Real, size: usize) -> Vec<Real> {
    // create the container for the data
    let mut data = Vec::with_capacity(size);

    // fill the container with the data
    for i in 0..size {
        data.push(start + (end - start) / (size - 1) as Real * i as Real);
    }

    return data; // return the data
}

pub fn repeat<T: NumberTrait>(vec: Vec<T>, n: usize) -> Vec<T> {
    (0..n).flat_map(|_| vec.iter().cloned()).collect()
}

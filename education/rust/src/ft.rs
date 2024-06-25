use complex::Complex;
use tensor::Tensor;
use number::Number;

pub fn dft<T: Number + From<f64> + Into<f64>>(input: &Tensor<Complex<T>>) -> Tensor<Complex<T>> {
    // define the output vector
    let mut output = Tensor::<Complex<T>>::init((&[input.len()]).to_vec());

    // compute the DFT
    for i in 0..input.len() {
        for j in 0..input.len() {
            output[i] = output[i] + input[j] * (-2.0 * Complex::i() * std::f64::consts::PI * i as f64 * j as f64 / input.len() as f64).exp();
        }
    }

    output // return the output vector
}

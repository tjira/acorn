mod fourier; mod number; mod vector;

use number::Real; use number::Complex; use vector::Vector;

fn main() {
    let a = Vector::new(vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0), Complex::new(4.0, 0.0)]);

    println!("{:?}", a);

    let a = fourier::fftcc(a);

    println!("{:?}", a);
}

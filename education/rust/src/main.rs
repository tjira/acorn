mod complex; mod ft; mod number; mod tensor;

use complex::Complex;
use tensor::Tensor;

fn main() {
    // let z = Complex::new(3.0, 4.0);

    // let v = Tensor::new(vec![Complex::new(1.0, 0.0), Complex::new(3.0, 0.0), Complex::new(5.0, 0.0), Complex::new(2.0, 0.0), Complex::new(9.0, 0.0)], vec![5]);
    let i = Tensor::linspace(Complex::new(0.0, 0.0), Complex::new(10.0, 0.0), vec![2, 3]);
    let fv = ft::dft(&i);

    // println!("{}", z);
    println!("{}", i);
    println!("{}", fv);

    // for (i, x) in fv.iter().enumerate() {
        // println!("f[{}] = {}", i, x);
    // }
}

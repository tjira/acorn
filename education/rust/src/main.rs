mod array; mod fourier; mod number; mod vector; mod wavefunction;

use number::Real; use number::Complex; use array::Array;

use vector::linspace; use vector::fftfreq;

use wavefunction::Wavefunction;

fn main() {
    // define the r-space and k-space arrays
    let r = Array::new(linspace(-8.0, 8.0, 64)); let ksq = Array::new(fftfreq(r.size(), r[1] - r[0])).apply(|x| x * x).cast::<Complex>();

    // define the potential array
    let u = r.clone().apply(|x| 0.5 * x * x).cast::<Complex>();

    // define the initial wavefunction
    let mut w = Wavefunction::new(&r, |x| (-x * x).exp(), 1.0).normalized();

    // define the propagators
    let (rp, kp) = w.propagators(&ksq, &u, 0.1);

    // propagation loop
    for _ in 0..1000 {

        // propagate the wavefunction
        w = w.propagated(&rp, &kp).normalized();

        // print the energy
        println!("{:?}", w.energy(&ksq, &u).norm());
    }
}

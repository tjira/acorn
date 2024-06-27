// crate includes
use num::complex::{Complex}; use ndarray::{Array, ArrayView, Axis, Ix1, Ix2, concatenate, prelude::array};

// personal includes
use fourier::{fftfreq}; use state::State; use wavefunction::Wavefunction;

pub struct QuantumDynamics {
    pub d: usize, pub iters: usize, pub nstate: usize, pub points: usize, pub rmin: f64, pub rmax: f64, pub imaginary: bool, pub savewfn: bool
}

impl QuantumDynamics {
    fn coords(&self) -> (Array::<f64, Ix2>, Array::<Complex<f64>, Ix1>) {
        // define the r-space and k-space array in one dimension
        let rr = Array::<f64, Ix1>::linspace(self.rmin, self.rmax, self.points); let kk = fftfreq(rr.len(), rr[1] - rr[0]);

        // define the r-space and k-space array in d dimensions
        let mut r = Array::<f64, Ix2>::zeros((self.points.pow(self.d as u32), self.d)); let mut k = Array::<f64, Ix2>::zeros((self.points.pow(self.d as u32), self.d));

        // fill the r-space and k-space arrays
        for i in 0..r.shape()[0] {
            for j in 0..r.shape()[1] {
                r[[i, self.d - j - 1]] = rr[i / self.points.pow(j as u32) % self.points];
                k[[i, self.d - j - 1]] = kk[i / self.points.pow(j as u32) % self.points];
            }
        }

        // define the k-space array squared and return
        let ksq = k.axis_iter(Axis(0)).map(|x| Complex::from(x.mapv(|x| x * x).sum())).collect::<Array::<Complex<f64>, Ix1>>(); (r, ksq)
    }
    pub fn run<U: Fn(ArrayView<f64, Ix1>) -> f64 + Clone, W: Fn(ArrayView<f64, Ix1>) -> f64 + Clone>(&self, ufv: Vec<U>, wfv: Vec<W>) -> (Array<f64, Ix2>, Option<Vec<Array<f64, Ix2>>>, Vec<f64>) {
        // define the r-space and k-space squared and the potential array
        let (r, ksq) = self.coords(); let u = r.axis_iter(Axis(0)).map(ufv[0].clone()).collect::<Array<f64, Ix1>>();

        // define the initial wavefunction for every state
        let mut states = vec![State::new(&r, wfv[0].clone(), 1.0).normalized(); self.nstate];

        // define the wavefunction and energy containers
        let mut wt = vec![array![[]]; self.nstate]; let mut ei = vec![0.0; self.nstate];

        // loop over each state
        for i in 0..self.nstate {

            // save the initial wavefunction if the flag provided
            if self.savewfn {wt[i] = concatenate(Axis(1), &[r.clone().view(), states[i].matrix().view()]).unwrap()}

            // define the propagators
            let (rp, kp) = states[i].propagators(&ksq, &u, if self.imaginary {Complex::<f64>::from(1.0)} else {Complex::<f64>::i()}, 0.1);

            // print the header
            println!("{:>6} {:>20}", "ITER", "ENERGY");

            // state propagation loop
            for j in 0..self.iters {

                // propagate the wavefunction and normalize if requested
                states[i] = states[i].clone().propagated(&rp, &kp); if self.imaginary {states[i] = states[i].clone().normalized()}

                // subtract lower eigenstates
                if self.imaginary {for k in 0..i {states[i] = states[i].clone() - states[k].clone() * states[k].overlap(&states[i])}}

                // save the wavefunction to the container
                if self.savewfn {wt[i] = concatenate(Axis(1), &[wt[i].view(), states[i].matrix().view()]).unwrap()}

                // print the iteration info
                println!("{:>6} {:>20.14}", j + 1, states[i].hamilton(&ksq, &u, &states[i]).re)
            }

            // assign energy and print empty line
            ei[i] = Wavefunction::new(vec![states[i].clone()]).hamilton(&ksq, &u).re; println!("")
        }

        // return the data
        (concatenate(Axis(1), &[r.view(), u.clone().into_shape((u.len(), 1)).unwrap().view()]).unwrap(), if self.savewfn {Some(wt)} else {None}, ei)
    }
}

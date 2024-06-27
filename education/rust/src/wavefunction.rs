// crate includes
use num::complex::Complex; use ndarray::{Array, Zip, Ix1};

// personal includes
use state::State;

#[derive(Clone, Debug)]
pub struct Wavefunction {
    states: Vec<State>
}

impl Wavefunction {
    // constructors for the wavefunction
    pub fn new(states: Vec<State>) -> Wavefunction {
        Wavefunction{states: states}
    }

    // wavefunction operators and functionals
    pub fn hamilton(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<f64, Ix1>) -> Complex<f64> {
        Zip::from(&self.states).and(&self.states).map_collect(|x, y| x.hamilton(ksq, u, &y)).sum()
    }
}

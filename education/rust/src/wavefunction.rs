use std::ops; use num::complex::Complex; use ndarray::{Array, ArrayView, Axis, Zip, Ix1, Ix2}; use fourier::{fftcc, ifftcc};

#[derive(Clone, Debug)]
pub struct State {
    data: Array<Complex<f64>, Ix1>, dr: f64, m: f64, d: usize
}

#[derive(Clone, Debug)]
pub struct Wavefunction {
    states: Vec<State>
}

impl State {
    // constructors for the state
    pub fn new<F: Fn(ArrayView<f64, Ix1>) -> f64>(r: &Array<f64, Ix2>, f: F, m: f64) -> State {
        State{data: r.axis_iter(Axis(0)).map(f).collect::<Array<f64, Ix1>>().mapv(|x| Complex::from(x)), dr: r[[1, r.shape()[1] - 1]] - r[[0, r.shape()[1] - 1]], m: m, d: r.shape()[1]}
    }

    // state operators and functionals
    pub fn hamilton(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<Complex<f64>, Ix1>, ket: &State) -> Complex<f64> {
        self.kin(ksq, &ket) + self.pot(u, &ket)
    }
    pub fn kin(&self, ksq: &Array<Complex<f64>, Ix1>, ket: &State) -> Complex<f64> {
        self.data.mapv(|x| x.conj()).dot(&ifftcc(ksq * fftcc(ket.data.clone(), self.d), self.d)) * self.dr / (2.0 * self.m)
    }
    pub fn pot(&self, u: &Array<Complex<f64>, Ix1>, ket: &State) -> Complex<f64> {
        Zip::from(&self.data).and(&ket.data).map_collect(|x, y| x.conj() * y).dot(u) * self.dr
    }
    pub fn normalized(self) -> State {
        let norm = self.overlap(&self).norm().sqrt(); State{data: self.data.mapv(|x| x / norm), dr: self.dr, m: self.m, d: self.d}
    }
    pub fn overlap(&self, ket: &State) -> Complex<f64> {
        Zip::from(&self.data).and(&ket.data).map_collect(|x, y| x.conj() * y).sum() * self.dr
    }
    pub fn propagated(self, rp: &Array<Complex<f64>, Ix1>, kp: &Array<Complex<f64>, Ix1>) -> State {
        State{data: rp * ifftcc(kp * fftcc(rp * self.data, self.d), self.d), dr: self.dr, m: self.m, d: self.d}
    }
    pub fn propagators(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<Complex<f64>, Ix1>, dt: f64) -> (Array<Complex<f64>, Ix1>, Array<Complex<f64>, Ix1>) {
        (u.mapv(|x| (-0.5 * x * dt).exp()), ksq.mapv(|x| (-0.5 * x * dt / self.m).exp()))
    }
}

impl Wavefunction {
    // constructors for the wavefunction
    pub fn new(states: Vec<State>) -> Wavefunction {
        Wavefunction{states: states}
    }

    // wavefunction operators and functionals
    pub fn hamilton(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<Complex<f64>, Ix1>) -> Complex<f64> {
        Zip::from(&self.states).and(&self.states).map_collect(|x, y| x.hamilton(ksq, u, &y)).sum()
    }
}

impl ops::Add<State> for State {
    type Output = State;
    fn add(self, rhs: State) -> State {
        State{data: self.data + rhs.data, dr: self.dr, m: self.m, d: self.d}
    }
}

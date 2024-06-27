use num::complex::Complex; use ndarray::{Array, ArrayView, Axis, Ix1, Ix2}; use fourier::{fftcc, ifftcc};

#[derive(Clone, Debug)]
pub struct Wavefunction {
    data: Array<Complex<f64>, Ix1>, dr: f64, m: f64, d: usize
}

impl Wavefunction {
    // constructors for the wavefunction
    pub fn new<F: Fn(ArrayView<f64, Ix1>) -> f64>(r: &Array<f64, Ix2>, f: F, m: f64) -> Wavefunction {
        Wavefunction{data: r.axis_iter(Axis(0)).map(f).collect::<Array<f64, Ix1>>().mapv(|x| Complex::from(x)), dr: r[[1, r.shape()[1] - 1]] - r[[0, r.shape()[1] - 1]], m: m, d: r.shape()[1]}
    }

    // wavefunction operators and functionals
    pub fn energy(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<Complex<f64>, Ix1>) -> Complex<f64> {
        self.ekin(ksq) + self.epot(u)
    }
    pub fn ekin(&self, ksq: &Array<Complex<f64>, Ix1>) -> Complex<f64> {
        self.data.mapv(|x| x.conj()).dot(&ifftcc(ksq * fftcc(self.data.clone(), self.d), self.d)) * self.dr / (2.0 * self.m)
    }
    pub fn epot(&self, u: &Array<Complex<f64>, Ix1>) -> Complex<f64> {
        self.data.mapv(|x| x * x.conj()).dot(u) * self.dr
    }
    pub fn normalized(self) -> Wavefunction {
        let norm = self.norm(); Wavefunction{data: self.data.mapv(|x| x / norm), dr: self.dr, m: self.m, d: self.d}
    }
    pub fn norm(&self) -> f64 {
        self.overlap().norm().sqrt()
    }
    pub fn overlap(&self) -> Complex<f64> {
        self.data.mapv(|x| x * x.conj()).sum() * self.dr
    }
    pub fn propagated(&self, rp: &Array<Complex<f64>, Ix1>, kp: &Array<Complex<f64>, Ix1>) -> Wavefunction {
        Wavefunction{data: rp * ifftcc(kp * fftcc(rp * &self.data, self.d), self.d), dr: self.dr, m: self.m, d: self.d}
    }
    pub fn propagators(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<Complex<f64>, Ix1>, dt: f64) -> (Array<Complex<f64>, Ix1>, Array<Complex<f64>, Ix1>) {
        (u.mapv(|x| (-0.5 * x * dt).exp()), ksq.mapv(|x| (-0.5 * x * dt / self.m).exp()))
    }
}

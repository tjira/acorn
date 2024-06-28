// system includes
use std::ops;

// crate includes
use num::complex::Complex; use ndarray::{Array, ArrayView, Axis, Zip, Ix1, Ix2, concatenate};

// personal includes
use fourier::{fftcc, ifftcc};

#[derive(Clone, Debug)]
pub struct Wavefunction {
    data: Array<Complex<f64>, Ix1>, dr: f64, m: f64, d: usize
}

impl Wavefunction {
    // constructors for the state
    pub fn new<F: Fn(ArrayView<f64, Ix1>) -> f64>(r: &Array<f64, Ix2>, f: F, m: f64) -> Wavefunction {
        Wavefunction{data: r.axis_iter(Axis(0)).map(f).collect::<Array<f64, Ix1>>().mapv(|x| Complex::from(x)), dr: r[[1, r.shape()[1] - 1]] - r[[0, r.shape()[1] - 1]], m: m, d: r.shape()[1]}
    }

    // getters
    pub fn matrix(&self) -> Array<f64, Ix2> {
        concatenate(Axis(1), &[self.data.mapv(|x| x.re).into_shape((self.data.len(), 1)).unwrap().view(), self.data.mapv(|x| x.im).into_shape((self.data.len(), 1)).unwrap().view()]).unwrap()
    }

    // state operators and functionals
    pub fn hamilton(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<f64, Ix1>, ket: &Wavefunction) -> Complex<f64> {
        self.kin(ksq, &ket) + self.pot(u, &ket)
    }
    pub fn kin(&self, ksq: &Array<Complex<f64>, Ix1>, ket: &Wavefunction) -> Complex<f64> {
        self.data.mapv(|x| x.conj()).dot(&ifftcc(ksq * fftcc(ket.data.clone(), self.d), self.d)) * self.dr / (2.0 * self.m)
    }
    pub fn pot(&self, u: &Array<f64, Ix1>, ket: &Wavefunction) -> Complex<f64> {
        (Zip::from(&self.data).and(&ket.data).map_collect(|x, y| x.conj() * y) * u).sum() * self.dr
    }
    pub fn normalized(self) -> Wavefunction {
        let norm = self.overlap(&self).norm().sqrt(); Wavefunction{data: self.data.mapv(|x| x / norm), dr: self.dr, m: self.m, d: self.d}
    }
    pub fn overlap(&self, ket: &Wavefunction) -> Complex<f64> {
        Zip::from(&self.data).and(&ket.data).map_collect(|x, y| x.conj() * y).sum() * self.dr
    }
    pub fn propagated(self, rp: &Array<Complex<f64>, Ix1>, kp: &Array<Complex<f64>, Ix1>) -> Wavefunction {
        Wavefunction{data: rp * ifftcc(kp * fftcc(rp * self.data, self.d), self.d), dr: self.dr, m: self.m, d: self.d}
    }
    pub fn propagators(&self, ksq: &Array<Complex<f64>, Ix1>, u: &Array<f64, Ix1>, unit: Complex<f64>, dt: f64) -> (Array<Complex<f64>, Ix1>, Array<Complex<f64>, Ix1>) {
        (u.mapv(|x| (-0.5 * x * unit * dt).exp()), ksq.mapv(|x| (-0.5 * x * unit * dt / self.m).exp()))
    }
}

impl ops::Add<Wavefunction> for Wavefunction {
    type Output = Wavefunction;
    fn add(self, rhs: Wavefunction) -> Wavefunction {
        Wavefunction{data: self.data + rhs.data, dr: self.dr, m: self.m, d: self.d}
    }
}

impl ops::Sub<Wavefunction> for Wavefunction {
    type Output = Wavefunction;
    fn sub(self, rhs: Wavefunction) -> Wavefunction {
        Wavefunction{data: self.data - rhs.data, dr: self.dr, m: self.m, d: self.d}
    }
}

impl ops::Mul<Complex<f64>> for Wavefunction {
    type Output = Wavefunction;
    fn mul(self, rhs: Complex<f64>) -> Wavefunction {
        Wavefunction{data: self.data * rhs, dr: self.dr, m: self.m, d: self.d}
    }
}

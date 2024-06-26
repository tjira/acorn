use number::Real; use number::Complex; use array::Array; use fourier::{fftcc, ifftcc};

pub struct Wavefunction {
    pub data: Array<Complex>, pub dr: Real, pub m: Real
}

impl Wavefunction {
    // constructors for the wavefunction
    pub fn new<F: Fn(Real) -> Real>(r: &Array<Real>, f: F, m: Real) -> Wavefunction {
        Wavefunction{data: r.clone().apply(f).cast::<Complex>(), dr: r[1] - r[0], m: m}
    }

    // wavefunction operators and functionals
    pub fn energy(&self, ksq: &Array<Complex>, u: &Array<Complex>) -> Complex {
        self.ekin(ksq) + self.epot(u)
    }
    pub fn ekin(&self, ksq: &Array<Complex>) -> Complex {
        self.data.clone().conj().dot(&ifftcc(ksq.clone() * fftcc(self.data.clone()))) * self.dr / (2.0 * self.m)
    }
    pub fn epot(&self, u: &Array<Complex>) -> Complex {
        self.data.clone().apply(|x| x * x.conj()).dot(u) * self.dr
    }
    pub fn mass(&self) -> Real {
        self.m
    }
    pub fn normalized(self) -> Wavefunction {
        let norm = self.norm(); Wavefunction{data: self.data.apply(|x| x / norm), dr: self.dr, m: self.m}
    }
    pub fn norm(&self) -> Real {
        self.overlap().norm().sqrt()
    }
    pub fn overlap(&self) -> Complex {
        self.data.clone().apply(|x| x * x.conj()).sum() * self.dr
    }
    pub fn propagated(&self, rp: &Array<Complex>, kp: &Array<Complex>) -> Wavefunction {
        Wavefunction{data: rp.clone() * ifftcc(kp.clone() * fftcc(rp.clone() * self.data.clone())), dr: self.dr, m: self.m}
    }
    pub fn propagators(&self, ksq: &Array<Complex>, u: &Array<Complex>, dt: Real) -> (Array<Complex>, Array<Complex>) {
        (u.clone().apply(|x| (-0.5 * x * dt).exp()), ksq.clone().apply(|x| (-0.5 * x * dt / self.m).exp()))
    }
}

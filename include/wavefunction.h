#pragma once

#include "linalg.h"

class Wavefunction {
public:
    Wavefunction(const std::string& guess, const EigenMatrix<>& r, double mass, double momentum = 0);

    // operators with other wavefunctions
    Wavefunction operator-(const Wavefunction& wfn) const {Wavefunction res(*this); res.data -= wfn.data; return res;}

    // operators with real scalars
    Wavefunction operator/(double scalar) const {Wavefunction res(*this); res.data /= scalar; return res;}

    // operators with complex scalars
    Wavefunction operator*(const std::complex<double>& scalar) const {Wavefunction res(*this); res.data *= scalar; return res;}

    // propagators and  other quantum operators
    double energy(const EigenMatrix<>& U) const; double norm() const; std::complex<double> overlap(const Wavefunction& wfn) const;
    Wavefunction propagate(const EigenMatrix<std::complex<double>>& R, const EigenMatrix<std::complex<double>>& K) const;

    // getters and input/output functions
    EigenMatrix<> getk() const {return k;} EigenMatrix<> getr() const {return r;} void Write(const std::string& path) const; EigenMatrix<std::complex<double>>& get() {return data;}

private:
    EigenMatrix<std::complex<double>> data; EigenMatrix<> r, k; double mass, dr;
};

#pragma once

#include "linalg.h"

template <int S> using MatrixOfMatrices = std::array<std::array<ComplexMatrix, S>, S>;

template <int S>
class Wavefunction {
public:
    Wavefunction(const Matrix& data, const Matrix& r, double mass = 1, double momentum = 0); Wavefunction() = default;

    // arithmetic operators
    Wavefunction operator*(const std::complex<double>& scalar) const {Wavefunction wfn(*this); for (int i = 0; i < S; i++) wfn.data.col(i) *= scalar; return wfn;}
    Wavefunction operator-(const Wavefunction& wfn) const {Wavefunction wfnn(*this); for (int i = 0; i < S; i++) wfnn.data.col(i) -= wfn.data.col(i); return wfnn;}

    // propagators and other quantum operators
    std::tuple<MatrixOfMatrices<S>, MatrixOfMatrices<S>> propagator(const Matrix& U, const std::complex<double>& unit, double step) const;
    double energy(const Matrix& U) const; std::complex<double> overlap(const Wavefunction& wfn) const; Vector norm() const;
    Wavefunction<S> propagate(const MatrixOfMatrices<S>& R, const MatrixOfMatrices<S>& K) const; Wavefunction normalized() const;
    Matrix density() const; Wavefunction<S> adiabatize(const std::vector<Matrix>& UT) const;

    // getters and input/output functions
    const ComplexMatrix& get() const {return data;} const Matrix& getr() const {return r;}

private:
    ComplexMatrix data; Matrix r, k; double mass, dr;
};

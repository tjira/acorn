#pragma once

#include "linalg.h"

template <int S> using MatrixOfMatrices = std::array<std::array<EigenMatrix<std::complex<double>>, S>, S>;

template <int S>
class Wavefunction {
public:
    Wavefunction(const EigenMatrix<>& data, const EigenMatrix<>& r, double mass = 1, double momentum = 0); Wavefunction() = default;

    // arithmetic operators
    Wavefunction operator*(const std::complex<double>& scalar) const {Wavefunction wfn(*this); for (int i = 0; i < S; i++) wfn.data.col(i) *= scalar; return wfn;}
    Wavefunction operator-(const Wavefunction& wfn) const {Wavefunction wfnn(*this); for (int i = 0; i < S; i++) wfnn.data.col(i) -= wfn.data.col(i); return wfnn;}

    // propagators and other quantum operators
    std::tuple<MatrixOfMatrices<S>, MatrixOfMatrices<S>> propagator(const EigenMatrix<>& U, const std::complex<double>& unit, double step) const;
    double energy(const EigenMatrix<>& U) const; std::complex<double> overlap(const Wavefunction& wfn) const; EigenVector<> norm() const;
    Wavefunction<S> propagate(const MatrixOfMatrices<S>& R, const MatrixOfMatrices<S>& K) const; Wavefunction normalized() const;
    EigenMatrix<> density() const; Wavefunction<S> adiabatize(const std::vector<EigenMatrix<>>& UT) const;

    // getters and input/output functions
    const EigenMatrix<std::complex<double>>& get() const {return data;} const EigenMatrix<>& getr() const {return r;}

private:
    EigenMatrix<std::complex<double>> data; EigenMatrix<> r, k; double mass, dr;
};

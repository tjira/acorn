#pragma once

#include "linalg.h"

namespace Acorn {
    namespace QDYN {
        class Wavefunction {
        public:
            Wavefunction(const Matrix& data, const Matrix& r, double mass = 1, double momentum = 0); Wavefunction();

            // arithmetic operators
            Wavefunction operator-(const Wavefunction& other) const; Wavefunction operator*(const std::complex<double>& scalar) const;

            // propagators and other quantum operators
            std::tuple<std::vector<ComplexMatrix>, std::vector<ComplexMatrix>> propagator(const Matrix& U, const std::complex<double>& unit, double step) const;
            Wavefunction propagate(const std::vector<ComplexMatrix>& R, const std::vector<ComplexMatrix>& K) const; Wavefunction normalized() const;
            double energy(const Matrix& U) const; std::complex<double> overlap(const Wavefunction& wfn) const;
            Matrix density() const; Wavefunction adiabatize(const std::vector<Matrix>& UT) const;

            // getters
            const ComplexMatrix& get() const; const Matrix& getr() const;

        private:
            ComplexMatrix data; Matrix r, k; std::vector<int> shape; double mass, dr;
        };
    }
}

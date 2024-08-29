#pragma once

#include "fourier.h"

namespace Acorn {
    namespace QDYN {
        class Wavefunction {
        public:
            Wavefunction(const Eigen::MatrixXd& data, const Eigen::MatrixXd& r, double mass = 1, double momentum = 0); Wavefunction();

            // arithmetic operators
            Wavefunction operator-(const Wavefunction& other) const; Wavefunction operator*(const std::complex<double>& scalar) const;

            // propagators and other quantum operators
            std::tuple<std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>> propagator(const Eigen::MatrixXd& U, const std::complex<double>& unit, double step) const;
            Wavefunction propagate(const std::vector<Eigen::MatrixXcd>& R, const std::vector<Eigen::MatrixXcd>& K) const; Wavefunction normalized() const;
            double energy(const Eigen::MatrixXd& U) const; std::complex<double> overlap(const Wavefunction& wfn) const;
            Eigen::MatrixXd density() const; Wavefunction adiabatize(const std::vector<Eigen::MatrixXd>& UT) const;

            // getters
            const Eigen::MatrixXcd& get() const; const Eigen::MatrixXd& getr() const;

        private:
            Eigen::MatrixXcd data; Eigen::MatrixXd r, k; std::vector<int> shape; double mass, dr;
        };
    }
}

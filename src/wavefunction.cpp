#include "fourier.h"
#include "wavefunction.h"

template <int D> const ComplexMatrix& Wavefunction<D>::get() const {return data;} template <int D> int Wavefunction<D>::states() const {return data.cols();} template <int D> const Matrix& Wavefunction<D>::getr() const {return r;}

template <int D>
Wavefunction<D> Wavefunction<D>::operator*(const std::complex<double>& scalar) const {
    Wavefunction<D> wfn(*this); wfn.data *= scalar; return wfn;
}

template <int D>
Wavefunction<D> Wavefunction<D>::operator-(const Wavefunction<D>& other) const {
    Wavefunction<D> wfn(*this); wfn.data -= other.data; return wfn;
}

template <int D>
Wavefunction<D>::Wavefunction(const Matrix& data, const Matrix& r, double mass, double momentum) : data(data.rows(), data.cols() / 2), r(r), k(r.rows(), r.cols()), mass(mass) {
    // add the momentum to the wavefunction
    for (int i = 0; i < 2 * this->data.cols(); i += 2) {
        this->data.col(i / 2) = (data.col(i) + std::complex<double>(0, 1) * data.col(i + 1)).array() * (std::complex<double>(0, 1) * momentum * r.array()).exp();
    } dr = r(1, r.cols() - 1) - r(0, r.cols() - 1);

    // calculate the k-space grid
    k.fill(2 * M_PI / k.size() / dr); for (int i = 0; i < k.size(); i++) k(i) *= i - (i < k.size() / 2 ? 0 : k.size());
}

template <int D>
Wavefunction<D> Wavefunction<D>::adiabatize(const std::vector<Matrix>& UT) const {
    // define adiabatic wavefunction
    Wavefunction<D> wfn(*this);

    // loop over all points and transform the wavefunction
    for (int j = 0; j < data.rows(); j++) wfn.data.row(j) = UT.at(j).adjoint() * data.row(j).transpose();

    // return the wfn
    return wfn;
}

template <int D>
Matrix Wavefunction<D>::density() const {
    // return the density matrix
    return (data.transpose() * data.conjugate()).array().abs() * dr;
}

template <int D>
double Wavefunction<D>::energy(const Matrix& U) const {
    // define the kinetic and potential energy
    double Ek = 0, Ep = 0;

    // loop over all wavefunctions
    for (int i = 0, n = data.cols(); i < data.cols(); i++) {

        // calculate the kinetic and potential energy
        Ek += (data.col(i).adjoint() * FourierTransform::IFFT(k.array().pow(2) * FourierTransform::FFT(data.col(i)).array()))(0).real() * dr;
        for (int j = 0; j < n; j++) Ep += (data.col(i).adjoint() * (U.col(i * n + j).array() * data.col(j).array()).matrix())(0).real() * dr;
    }

    // return the total energy
    return 0.5 * Ek / mass + Ep;
}

template <int D>
Wavefunction<D> Wavefunction<D>::normalized() const {
    // copy the wavefunction and normalize it, then return it
    Wavefunction<D> wfn(*this); wfn.data /= std::sqrt(std::abs(overlap(wfn))); return wfn;
}

template <int D>
std::complex<double> Wavefunction<D>::overlap(const Wavefunction<D>& wfn) const {
    // define the overlap container
    std::complex<double> overlap = 0;
    
    // calculate all the contributions to the overlap
    for (int i = 0; i < data.cols(); i++) for (int j = 0; j < wfn.data.cols(); j++) overlap += (data.col(i).adjoint() * wfn.data.col(j))(0) * dr;

    // return
    return overlap;
}

template <int D>
Wavefunction<D> Wavefunction<D>::propagate(const std::vector<ComplexMatrix>& R, const std::vector<ComplexMatrix>& K) const {
    // output wavefunction
    Wavefunction wfn(*this);

    // perform the half step in the real space
    for (int i = 0; i < data.rows(); i++) wfn.data.row(i) = R.at(i) * wfn.data.row(i).transpose();

    // perform the fourier transform
    for (int j = 0; j < data.cols(); j++) wfn.data.col(j) = FourierTransform::FFT(wfn.data.col(j));

    // perform the full step in the momentum space
    for (int i = 0; i < data.rows(); i++) wfn.data.row(i) = K.at(i) * wfn.data.row(i).transpose();

    // perform the inverse fourier transform
    for (int j = 0; j < data.cols(); j++) wfn.data.col(j) = FourierTransform::IFFT(wfn.data.col(j));

    // perform the half step in the real space
    for (int i = 0; i < data.rows(); i++) wfn.data.row(i) = R.at(i) * wfn.data.row(i).transpose();

    // return
    return wfn;
}

template <int D>
std::tuple<std::vector<ComplexMatrix>, std::vector<ComplexMatrix>> Wavefunction<D>::propagator(const Matrix& U, const std::complex<double>& unit, double step) const {
    // define the propagator array
    std::vector<ComplexMatrix> R(U.rows()), K(U.rows());

    // calculate the propagators
    for (int i = 0, n = data.cols(); i < U.rows(); i++) {

        // initialize the eigenvalue solver and diagonalize the potential
        Eigen::SelfAdjointEigenSolver<Matrix> solver(U.row(i).reshaped(n, n));
        ComplexMatrix C = solver.eigenvectors(), E = solver.eigenvalues();

        // create the real and momentum space propagator matrix for the current point
        ComplexMatrix KI = (-0.5 * unit * std::pow(k(i), 2) * step * Vector::Ones(n).array() / mass).exp();
        ComplexMatrix RI = C * (-0.5 * unit * E.array() * step).exp().matrix().asDiagonal() * C.adjoint();

        // store the propagator matrices
        R.at(i) = RI, K.at(i) = KI.asDiagonal();
    }

    // return the propagator
    return {R, K};
}

template class Wavefunction<1>;

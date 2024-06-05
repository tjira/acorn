#include "fourier.h"
#include "wavefunction.h"

template<int S>
Wavefunction<S>::Wavefunction(const Matrix& data, const Matrix& r, double mass, double momentum) : data(data.rows(), data.cols() / 2), r(r), k(r.rows(), r.cols()), mass(mass) {
    // add the momentum to the wavefunction
    for (int i = 0; i < 2 * S; i += 2) {
        this->data.col(i / 2) = (data.col(i) + std::complex<double>(0, 1) * data.col(i + 1)).array() * (std::complex<double>(0, 1) * momentum * r.array()).exp();
    } dr = r(1) - r(0);

    // calculate the k-space grid
    k.fill(2 * M_PI / k.size() / dr); for (int i = 0; i < k.size(); i++) k(i) *= i - (i < k.size() / 2 ? 0 : k.size());
}

template <int S>
Wavefunction<S> Wavefunction<S>::adiabatize(const std::vector<Matrix>& UT) const {
    // define adiabatic wavefunction
    Wavefunction wfn(*this);

    // loop over all points and transform the wavefunction
    for (int j = 0; j < data.rows(); j++) wfn.data.row(j) = UT.at(j).adjoint() * data.row(j).transpose();

    // return the wfn
    return wfn;
}

template <int S>
Matrix Wavefunction<S>::density() const {
    // return the density matrix
    return (data.transpose() * data.conjugate()).array().abs() * dr;
}

template <>
double Wavefunction<1>::energy(const Matrix& U) const {
    // calculate the kinetic energy
    double Ek = (data.adjoint() * FourierTransform::IFFT(k.array().pow(2) * FourierTransform::FFT(data).array()))(0).real() * dr;

    // calculate the potential energy
    double Ep = (data.adjoint() * (U.array() * data.array()).matrix())(0).real() * dr;

    // return the total energy
    return 0.5 * Ek / mass + Ep;
}

template <>
double Wavefunction<2>::energy(const Matrix& U) const {
    // calculate the kinetic energies
    double Ek00 = (data.col(0).adjoint() * FourierTransform::IFFT(k.array().pow(2) * FourierTransform::FFT(data.col(0)).array()))(0).real() * dr;
    double Ek11 = (data.col(1).adjoint() * FourierTransform::IFFT(k.array().pow(2) * FourierTransform::FFT(data.col(1)).array()))(0).real() * dr;

    // calculate the potential energies
    double Ep00 = (data.col(0).adjoint() * (U.col(0).array() * data.col(0).array()).matrix())(0).real() * dr;
    double Ep01 = (data.col(0).adjoint() * (U.col(1).array() * data.col(1).array()).matrix())(0).real() * dr;
    double Ep10 = (data.col(1).adjoint() * (U.col(2).array() * data.col(0).array()).matrix())(0).real() * dr;
    double Ep11 = (data.col(1).adjoint() * (U.col(3).array() * data.col(1).array()).matrix())(0).real() * dr;

    // return the total energy
    return 0.5 * (Ek00 + Ek11) / mass + Ep00 + Ep01 + Ep10 + Ep11;
}

template <int S>
Vector Wavefunction<S>::norm() const {
    // initialize the vector
    Vector norm(S);

    // calculate the norm of the wavefunction and store it in the vector
    for (int i = 0; i < S; i++) norm(i) = std::sqrt(data.col(i).array().abs2().sum() * dr);

    // return the norm
    return norm;
}

template <int S>
Wavefunction<S> Wavefunction<S>::normalized() const {
    // copy the wavefunction and calculate the norms
    Wavefunction wfn(*this); Vector norms(S);

    // calculate the norms of the wavefunction
    for (int i = 0; i < S; i++) norms(i) = std::sqrt(data.col(i).array().abs2().sum() * dr);

    // divide the wavefunction by the norm
    for (int i = 0; i < S; i++) if (norms(i) > 1e-14) wfn.data.col(i) /= norms(i);

    // return the wfn
    return wfn;
}

template <>
std::complex<double> Wavefunction<1>::overlap(const Wavefunction<1>& wfn) const {
    // return the overlap of the wavefunctions
    return (wfn.data.col(0).conjugate().array() * data.col(0).array()).sum() * dr;
}

template <>
Wavefunction<1> Wavefunction<1>::propagate(const MatrixOfMatrices<1>& R, const MatrixOfMatrices<1>& K) const {
    // define the next wfn
    Wavefunction wfn(*this);

    // propagate the wavefunction
    wfn.data.col(0) = R.at(0).at(0).array() * wfn.data.col(0).array();
    wfn.data.col(0) = FourierTransform::IFFT(K.at(0).at(0).array() * FourierTransform::FFT(wfn.data.col(0)).array());
    wfn.data.col(0) = R.at(0).at(0).array() * wfn.data.col(0).array();

    // return the wfn
    return wfn;
}

template <>
Wavefunction<2> Wavefunction<2>::propagate(const MatrixOfMatrices<2>& R, const MatrixOfMatrices<2>& K) const {
    // define the next and temp wfn
    Wavefunction wfn(*this), wfnt(*this);

    // propagate the wavefunction
    wfnt.data.col(0) = R.at(0).at(0).array() * wfn.data.col(0).array() + R.at(0).at(1).array() * wfn.data.col(1).array();
    wfnt.data.col(1) = R.at(1).at(0).array() * wfn.data.col(0).array() + R.at(1).at(1).array() * wfn.data.col(1).array(); wfn = wfnt;
    wfnt.data.col(0) = FourierTransform::IFFT(K.at(0).at(0).array() * FourierTransform::FFT(wfn.data.col(0)).array());
    wfnt.data.col(1) = FourierTransform::IFFT(K.at(1).at(1).array() * FourierTransform::FFT(wfn.data.col(1)).array()); wfn = wfnt;
    wfnt.data.col(0) = R.at(0).at(0).array() * wfn.data.col(0).array() + R.at(0).at(1).array() * wfn.data.col(1).array();
    wfnt.data.col(1) = R.at(1).at(0).array() * wfn.data.col(0).array() + R.at(1).at(1).array() * wfn.data.col(1).array(); wfn = wfnt;

    // return the wfn
    return wfn;
}

template <int S>
std::tuple<MatrixOfMatrices<S>, MatrixOfMatrices<S>> Wavefunction<S>::propagator(const Matrix& U, const std::complex<double>& unit, double step) const {
    // define the propagator array
    MatrixOfMatrices<S> R, K;

    // fill the momentum space propagator
    for (int i = 0; i < S; i++) K.at(i).at(i) = (-0.5 * unit * k.array().pow(2) * step / mass).exp();

    // initialize the real space propagator
    for (int i = 0; i < S; i++) for (int j = 0; j < S; j++) R.at(i).at(j) = Vector::Zero(U.rows());

    // calculate the real space propagator
    for (int i = 0; i < U.rows(); i++) {

        // initialize the eigenvalue solver and diagonalize the potential
        Eigen::SelfAdjointEigenSolver<Matrix> solver(U.row(i).reshaped(S, S));
        ComplexMatrix C = solver.eigenvectors(), D = solver.eigenvalues();

        // create the real space propagator matrix for the current point
        ComplexMatrix RI = C * (-0.5 * unit * D.array() * step).exp().matrix().asDiagonal() * C.adjoint();

        // store the propagator matrix
        for (int j = 0; j < S; j++) for (int k = 0; k < S; k++) R.at(j).at(k)(i) = RI(j, k);
    }

    // return the propagator
    return {R, K};
}

template class Wavefunction<1>;
template class Wavefunction<2>;

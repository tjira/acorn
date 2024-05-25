#include "fourier.h"
#include "wavefunction.h"

#include <iostream>

template<int S>
Wavefunction<S>::Wavefunction(const EigenMatrix<std::complex<double>>& data, const EigenMatrix<>& r, double mass, double momentum) : r(r), k(r.rows(), r.cols()), mass(mass) {
    // add the momentum to the wavefunction
    for (int i = 0; i < 2 * S; i += 2) {
        this->data.at(i / 2) = (data.col(i) + std::complex<double>(0, 1) * data.col(i + 1)).array() * (std::complex<double>(0, 1) * momentum * r.array()).exp();
    } dr = r(1) - r(0);

    // calculate the k-space grid
    k.fill(2 * M_PI / k.size() / dr); for (int i = 0; i < k.size(); i++) k(i) *= i - (i < k.size() / 2 ? 0 : k.size());
}

template <int S>
EigenMatrix<> Wavefunction<S>::density() const {
    // define the density matrix
    EigenMatrix<> density(S, S); density.setZero();

    for (int i = 0; i < S; i++) {
        for (int j = 0; j < S; j++) {
            density(i, j) += std::abs((data.at(i).transpose() * Eigen::ComplexConjugate(data.at(j)))(0)) * dr;
        }
    }

    // return the density matrix
    return density;
}

template <>
double Wavefunction<1>::energy(const EigenMatrix<>& U) const {
    // obtain the negative spatial second derivative of the wavefunction in the k-space
    EigenMatrix<std::complex<double>> wfnk = k.array().pow(2) * FourierTransform::Forward(data.at(0)).array();

    // calculate the kinetic energy
    double Ek = 0.5 * (Eigen::ComplexConjugate(data.at(0)).transpose() * FourierTransform::Inverse(wfnk))(0).real() * dr / mass;

    // calculate the potential energy
    double Ep = (Eigen::ComplexConjugate(data.at(0)).transpose() * (U.array() * data.at(0).array()).matrix())(0).real() * dr;

    // return the total energy
    return Ek + Ep;
}

template <>
double Wavefunction<2>::energy(const EigenMatrix<>& U) const {
    // calculate the kinetic energies
    EigenMatrix<std::complex<double>> Ek00 = Eigen::ComplexConjugate(data.at(0)).array() * FourierTransform::Inverse(k.array().pow(2) * FourierTransform::Forward(data.at(0)).array()).array();
    EigenMatrix<std::complex<double>> Ek11 = Eigen::ComplexConjugate(data.at(1)).array() * FourierTransform::Inverse(k.array().pow(2) * FourierTransform::Forward(data.at(1)).array()).array();

    // calculate the potential energies
    EigenMatrix<std::complex<double>> Ep00 = Eigen::ComplexConjugate(data.at(0)).array() * U.col(0).array() * data.at(0).array();
    EigenMatrix<std::complex<double>> Ep01 = Eigen::ComplexConjugate(data.at(0)).array() * U.col(1).array() * data.at(1).array();
    EigenMatrix<std::complex<double>> Ep10 = Eigen::ComplexConjugate(data.at(1)).array() * U.col(2).array() * data.at(0).array();
    EigenMatrix<std::complex<double>> Ep11 = Eigen::ComplexConjugate(data.at(1)).array() * U.col(3).array() * data.at(1).array();

    // return the total energy
    return (0.5 * (Ek00 + Ek11) / mass + Ep00 + Ep01 + Ep10 + Ep11).sum().real() * dr;
}

template <int S>
EigenVector<> Wavefunction<S>::norm() const {
    // initialize the vector
    EigenVector<> norm(S);

    // calculate the norm of the wavefunction and store it in the vector
    for (int i = 0; i < S; i++) norm(i) = std::sqrt(data.at(i).array().abs2().sum() * dr);

    // return the norm
    return norm;
}

template <int S>
Wavefunction<S> Wavefunction<S>::normalized() const {
    // copy the wavefunction and calculate the norms
    Wavefunction wfn(*this); EigenVector<> norms(S);

    // calculate the norms of the wavefunction
    for (int i = 0; i < S; i++) norms(i) = std::sqrt(data.at(i).array().abs2().sum() * dr);

    // divide the wavefunction by the norm
    for (int i = 0; i < S; i++) if (norms(i) > 1e-14) wfn.data.at(i) /= norms(i);

    // return the wfn
    return wfn;
}

template <>
std::complex<double> Wavefunction<1>::overlap(const Wavefunction<1>& wfn) const {
    // return the overlap of the wavefunctions
    return (Eigen::ComplexConjugate(wfn.data.at(0)).array() * data.at(0).array()).sum() * dr;
}

template <>
Wavefunction<1> Wavefunction<1>::propagate(const MatrixOfMatrices<1>& R, const MatrixOfMatrices<1>& K) const {
    // define the next wfn
    Wavefunction wfn(*this);

    // propagate the wavefunction
    wfn.data.at(0) = R.at(0).at(0).array() * wfn.data.at(0).array();
    wfn.data.at(0) = FourierTransform::Forward(wfn.data.at(0));
    wfn.data.at(0) = K.at(0).at(0).array() * wfn.data.at(0).array();
    wfn.data.at(0) = FourierTransform::Inverse(wfn.data.at(0));
    wfn.data.at(0) = R.at(0).at(0).array() * wfn.data.at(0).array();

    // return the wfn
    return wfn;
}

template <>
Wavefunction<2> Wavefunction<2>::propagate(const MatrixOfMatrices<2>& R, const MatrixOfMatrices<2>& K) const {
    // define the next and temp wfn
    Wavefunction wfn(*this), wfnt(*this);

    // propagate the wavefunction
    wfnt.data.at(0) = R.at(0).at(0).array() * wfn.data.at(0).array() + R.at(0).at(1).array() * wfn.data.at(1).array();
    wfnt.data.at(1) = R.at(1).at(0).array() * wfn.data.at(0).array() + R.at(1).at(1).array() * wfn.data.at(1).array();
    wfn = wfnt;
    wfnt.data.at(0) = FourierTransform::Inverse(K.at(0).at(0).array() * FourierTransform::Forward(wfn.data.at(0)).array());
    wfnt.data.at(1) = FourierTransform::Inverse(K.at(1).at(1).array() * FourierTransform::Forward(wfn.data.at(1)).array());
    wfn = wfnt;
    wfnt.data.at(0) = R.at(0).at(0).array() * wfn.data.at(0).array() + R.at(0).at(1).array() * wfn.data.at(1).array();
    wfnt.data.at(1) = R.at(1).at(0).array() * wfn.data.at(0).array() + R.at(1).at(1).array() * wfn.data.at(1).array();
    wfn = wfnt;

    // return the wfn
    return wfn;
}

template <>
std::tuple<MatrixOfMatrices<1>, MatrixOfMatrices<1>> Wavefunction<1>::propagator(const EigenMatrix<>& U, const std::complex<double>& unit, double step) const {
    // define the propagator array
    MatrixOfMatrices<1> R, K;

    // calculate the propagator
    R.at(0).at(0) = (-0.5 * unit * U.array() * step).exp(), K.at(0).at(0) = (-0.5 * unit * k.array().pow(2) * step / mass).exp();

    // return the propagator
    return {R, K};
}

template <>
std::tuple<MatrixOfMatrices<2>, MatrixOfMatrices<2>> Wavefunction<2>::propagator(const EigenMatrix<>& U, const std::complex<double>& unit, double step) const {
    // define the propagator array
    MatrixOfMatrices<2> R, K;

    // calculate the propagator
    K.at(0).at(0) = (-0.5 * unit * k.array().pow(2) * step / mass).exp();
    K.at(1).at(1) = (-0.5 * unit * k.array().pow(2) * step / mass).exp();
    EigenVector<std::complex<double>> D = 4 * U.col(2).array().abs().pow(2) + (U.col(0) - U.col(3)).array().pow(2);
    EigenVector<std::complex<double>> a = (-0.25 * unit * (U.col(0) + U.col(3)) * step).array().exp();
    EigenVector<std::complex<double>> b = (0.25 * D.array().sqrt() * step).cos();
    EigenVector<std::complex<double>> c = unit * (0.25 * D.array().sqrt() * step).sin() / D.array().sqrt();
    R.at(0).at(0) = a.array() * (b.array() + c.array() * (U.col(3) - U.col(0)).array());
    R.at(0).at(1) = -2 * a.array() * c.array() * U.col(1).array();
    R.at(1).at(0) = -2 * a.array() * c.array() * U.col(1).array();
    R.at(1).at(1) = a.array() * (b.array() + c.array() * (U.col(0) - U.col(3)).array());

    // return the propagator
    return {R, K};
}

template class Wavefunction<1>;
template class Wavefunction<2>;

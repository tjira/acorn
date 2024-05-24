#include "expression.h"
#include "fourier.h"
#include "wavefunction.h"

Wavefunction::Wavefunction(const Wavefunction& wfn) : data(wfn.data), r(wfn.r), k(wfn.k), mass(wfn.mass), dr(wfn.dr) {}

Wavefunction::Wavefunction(const std::string& guess, const EigenMatrix<>& r, double mass) : r(r), k(r.rows(), r.cols()), mass(mass), dr(r(1) - r(0)) {
    // evaluate the guess expression
    Expression wfnexpr(guess, {"x"}); data = wfnexpr.eval(r);

    // calculate the k-space grid
    k.fill(2 * M_PI / k.size() / dr); for (int i = 0; i < k.size(); i++) k(i) *= i - (i < k.size() / 2 ? 0 : k.size());
}

double Wavefunction::energy(const EigenMatrix<>& U) const {
    // obtain the negative spatial second derivative of the wavefunction in the k-space
    EigenMatrix<std::complex<double>> wfnk = k.array().pow(2) * FourierTransform::Forward(data).array();

    // calculate the kinetic energy
    double Ek = 0.5 * (Eigen::ComplexConjugate(data).transpose() * FourierTransform::Inverse(wfnk))(0).real() * dr / mass;

    // calculate the potential energy
    double Ep = (Eigen::ComplexConjugate(data).transpose() * (U.array() * data.array()).matrix())(0).real() * dr;

    // return the total energy
    return Ek + Ep;
}

double Wavefunction::norm() const {
    // return the norm of the wavefunction
    return std::sqrt(data.array().abs2().sum() * dr);
}

std::complex<double> Wavefunction::overlap(const Wavefunction& wfn) const {
    // return the overlap of the wavefunctions
    return (Eigen::ComplexConjugate(wfn.data).array() * data.array()).sum() * dr;
}

Wavefunction Wavefunction::propagate(const EigenMatrix<std::complex<double>>& R, const EigenMatrix<std::complex<double>>& K) const {
    // define the next wfn
    Wavefunction wfn(*this);

    // propagate the wavefunction
    wfn.data = R.array() * wfn.data.array();
    wfn.data = FourierTransform::Forward(wfn.data);
    wfn.data = K.array() * wfn.data.array();
    wfn.data = FourierTransform::Inverse(wfn.data);
    wfn.data = R.array() * wfn.data.array();

    // return the wfn
    return wfn;
}

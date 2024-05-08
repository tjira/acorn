#include "propagator.h"

template <int S>
Matrix<std::complex<double>> Propagator<S>::Propagate(const MatrixOfMatrices<Matrix<std::complex<double>>>& R, const MatrixOfMatrices<Matrix<std::complex<double>>>& K, Matrix<std::complex<double>> psi) {
    if constexpr (S == 1) {
        psi = R(0, 0).array() * psi.array();
        psi = Numpy::FFT(K(0, 0).array() * Numpy::FFT(psi).array(), 1);
        psi = R(0, 0).array() * psi.array();
    } else if constexpr (Matrix<std::complex<double>> psit = psi; S == 2) {
        psit.col(0) = R(0, 0).array() * psi.col(0).array() + R(0, 1).array() * psi.col(1).array();
        psit.col(1) = R(1, 0).array() * psi.col(0).array() + R(1, 1).array() * psi.col(1).array(); psi = psit;
        psit.col(0) = Numpy::FFT(K(0, 0).array() * Numpy::FFT(psi.col(0)).array(), 1);
        psit.col(1) = Numpy::FFT(K(1, 1).array() * Numpy::FFT(psi.col(1)).array(), 1); psi = psit;
        psit.col(0) = R(0, 0).array() * psi.col(0).array() + R(0, 1).array() * psi.col(1).array();
        psit.col(1) = R(1, 0).array() * psi.col(0).array() + R(1, 1).array() * psi.col(1).array(); psi = psit;
    } else {throw std::runtime_error("PROPAGATOR NOT IMPLEMENTED");} return psi;
}

template <int S>
double Propagator<S>::Energy(const ModelSystem& system, const MatrixOfMatrices<Matrix<>>& V, const Matrix<>& r, const Matrix<>& ksq, Matrix<std::complex<double>> psi) {
    if constexpr (S == 1) {
        Matrix<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * Numpy::FFT(ksq.array() * Numpy::FFT(psi).array(), 1).array();
        Matrix<std::complex<double>> Ep = EigenConj(psi).array() * V(0, 0).array() * psi.array();
        return (Ek / system.mass() + Ep).sum().real() * std::pow(r(1) - r(0), system.vars().size());
    } else if constexpr (S == 2) {
        Matrix<std::complex<double>> Ek00 = 0.5 * EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * Numpy::FFT(ksq.array() * Numpy::FFT(psi.col(0)).array(), 1).array();
        Matrix<std::complex<double>> Ek11 = 0.5 * EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * Numpy::FFT(ksq.array() * Numpy::FFT(psi.col(1)).array(), 1).array();
        Matrix<std::complex<double>> Ep00 = EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * V(0, 0).array() * psi.col(0).array();
        Matrix<std::complex<double>> Ep01 = EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * V(0, 1).array() * psi.col(1).array();
        Matrix<std::complex<double>> Ep10 = EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * V(1, 0).array() * psi.col(0).array();
        Matrix<std::complex<double>> Ep11 = EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * V(1, 1).array() * psi.col(1).array();
        return ((Ek00 + Ek11) / system.mass() + Ep00 + Ep01 + Ep10 + Ep11).sum().real() * std::pow(r(1) - r(0), system.vars().size());
    } else throw std::runtime_error("PROPAGATOR NOT IMPLEMENTED");
}

template <int S>
std::tuple<ComplexMatrixOfMatrices, ComplexMatrixOfMatrices> Propagator<S>::Get(const ModelSystem& system, const ComplexMatrixOfMatrices& V, const Matrix<>& ksq, double step, bool real) {
    std::complex<double> unit = real ? std::complex<double>(0, 1) : std::complex<double>(1, 0);
    if constexpr (MatrixOfMatrices<Matrix<std::complex<double>>> R(1, 1), K(1, 1); S == 1) {
        K(0, 0) = (-0.5 * unit * ksq.array() * step / system.mass()).exp();
        R(0, 0) = (-0.5 * unit * V(0, 0).array() * step).exp();
        return {R, K};
    } else if (MatrixOfMatrices<Matrix<std::complex<double>>> R(2, 2), K(2, 2); S == 2) {
        K(0, 0) = (-0.5 * unit * ksq.array() * step / system.mass()).exp();
        K(1, 1) = (-0.5 * unit * ksq.array() * step / system.mass()).exp();
        Vector<std::complex<double>> D = 4 * V(1, 0).array().abs().pow(2) + (V(0, 0) - V(1, 1)).array().pow(2);
        Vector<std::complex<double>> a = (-0.25 * unit * (V(0, 0) + V(1, 1)) * step).array().exp();
        Vector<std::complex<double>> b = (0.25 * D.array().sqrt() * step).cos();
        Vector<std::complex<double>> c = unit * (0.25 * D.array().sqrt() * step).sin() / D.array().sqrt();
        R(0, 0) = a.array() * (b.array() + c.array() * (V(1, 1) - V(0, 0)).array());
        R(0, 1) = -2 * a.array() * c.array() * V(0, 1).array();
        R(1, 0) = -2 * a.array() * c.array() * V(0, 1).array();
        R(1, 1) = a.array() * (b.array() + c.array() * (V(0, 0) - V(1, 1)).array());
        return {R, K};
    } else throw std::runtime_error("PROPAGATOR NOT IMPLEMENTED");
}

template class Propagator<1>;
template class Propagator<2>;

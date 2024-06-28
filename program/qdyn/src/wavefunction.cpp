#include "fourier.h"
#include "wavefunction.h"

Wavefunction::Wavefunction() = default; const ComplexMatrix& Wavefunction::get() const {return data;} const Matrix& Wavefunction::getr() const {return r;}

Wavefunction::Wavefunction(const Matrix& data, const Matrix& r, double mass, double momentum) : data(data.rows(), data.cols() / 2), r(r), k(r.rows(), r.cols()), mass(mass) {
    // add the momentum to the wavefunction
    for (int i = 0; i < 2 * this->data.cols(); i += 2) {
        this->data.col(i / 2) = (data.col(i) + std::complex<double>(0, 1) * data.col(i + 1)).array() * (std::complex<double>(0, 1) * momentum * r.rowwise().sum().array()).exp();
    } dr = r(1, r.cols() - 1) - r(0, r.cols() - 1);

    // get the grid size and create the shapes array, suppose the grid is square
    int n = std::round(std::pow(r.rows(), 1.0 / r.cols())); for (int i = 0; i < r.cols(); i++) shape.push_back(n);

    // calculate the last column of the k-space grid
    for (int i = 0; i < r.rows() / n; i++) {
        k.rightCols(1).topRows((i + 1) * n).bottomRows(n).fill(2 * M_PI / n / dr); for (int j = 0; j < n; j++) k.rightCols(1)(i * n + j) *= j - (j < n / 2 ? 0 : n);
    }

    // calculate the rest of the k-space grid
    for (int i = 0; i < r.cols() - 1; i++) for (int j = 0; j < r.rows(); j++) k(j, i) = k(j / (int)std::round(std::pow(n, r.cols() - i - 1)), r.cols() - 1);
}

Wavefunction Wavefunction::operator*(const std::complex<double>& scalar) const {
    Wavefunction wfn(*this); wfn.data *= scalar; return wfn;
}

Wavefunction Wavefunction::operator-(const Wavefunction& other) const {
    Wavefunction wfn(*this); wfn.data -= other.data; return wfn;
}

Wavefunction Wavefunction::adiabatize(const std::vector<Matrix>& UT) const {
    // define adiabatic wavefunction
    Wavefunction wfn(*this);

    // loop over all points and transform the wavefunction
    for (int j = 0; j < data.rows(); j++) wfn.data.row(j) = UT.at(j).adjoint() * data.row(j).transpose();

    // return the wfn
    return wfn;
}

Matrix Wavefunction::density() const {
    // return the density matrix
    return (data.transpose() * data.conjugate()).array().abs() * dr;
}

double Wavefunction::energy(const Matrix& U) const {
    // define the kinetic and potential energy
    double Ek = 0, Ep = 0;

    // loop over all wavefunctions
    for (int i = 0, n = data.cols(); i < data.cols(); i++) {

        // calculate the kinetic energy contribution
        Ek += (data.col(i).adjoint() * FourierTransform::IFFT(k.array().pow(2).rowwise().sum() * FourierTransform::FFT(data.col(i), shape).array(), shape))(0).real() * dr;

        // calculate the potential energy contributions
        for (int j = 0; j < n; j++) Ep += (data.col(i).adjoint() * (U.col(i * n + j).array() * data.col(j).array()).matrix())(0).real() * dr;
    }

    // return the total energy
    return 0.5 * Ek / mass + Ep;
}

Wavefunction Wavefunction::normalized() const {
    // copy the wavefunction and normalize it, then return it
    Wavefunction wfn(*this); wfn.data /= std::sqrt(std::abs(overlap(wfn))); return wfn;
}

std::complex<double> Wavefunction::overlap(const Wavefunction& wfn) const {
    // define the overlap container
    std::complex<double> overlap = 0;
    
    // calculate all the contributions to the overlap
    for (int i = 0; i < data.cols(); i++) for (int j = 0; j < wfn.data.cols(); j++) overlap += (data.col(i).adjoint() * wfn.data.col(j))(0) * dr;

    // return
    return overlap;
}

Wavefunction Wavefunction::propagate(const std::vector<ComplexMatrix>& R, const std::vector<ComplexMatrix>& K) const {
    // output wavefunction
    Wavefunction wfn(*this);

    // perform the half step in the real space
    for (int i = 0; i < data.rows(); i++) wfn.data.row(i) = R.at(i) * wfn.data.row(i).transpose();

    // perform the fourier transform
    for (int j = 0; j < data.cols(); j++) wfn.data.col(j) = FourierTransform::FFT(wfn.data.col(j), shape);

    // perform the full step in the momentum space
    for (int i = 0; i < data.rows(); i++) wfn.data.row(i) = K.at(i) * wfn.data.row(i).transpose();

    // perform the inverse fourier transform
    for (int j = 0; j < data.cols(); j++) wfn.data.col(j) = FourierTransform::IFFT(wfn.data.col(j), shape);

    // perform the half step in the real space
    for (int i = 0; i < data.rows(); i++) wfn.data.row(i) = R.at(i) * wfn.data.row(i).transpose();

    // return
    return wfn;
}

std::tuple<std::vector<ComplexMatrix>, std::vector<ComplexMatrix>> Wavefunction::propagator(const Matrix& U, const std::complex<double>& unit, double step) const {
    // define the propagator array
    std::vector<ComplexMatrix> R(U.rows()), K(U.rows());

    // calculate the propagators
    for (int i = 0, n = data.cols(); i < U.rows(); i++) {

        // initialize the eigenvalue solver and diagonalize the potential
        Eigen::SelfAdjointEigenSolver<Matrix> solver(U.row(i).reshaped(n, n));
        ComplexMatrix C = solver.eigenvectors(), E = solver.eigenvalues();

        // momentum space propagator as a diagonal matrix with kinetic energy operators
        ComplexMatrix KI = (-0.5 * unit * k.row(i).array().pow(2).sum() * step * Vector::Ones(n).array() / mass).exp();

        // real space propagator as a general matrix exponential
        ComplexMatrix RI = C * (-0.5 * unit * E.array() * step).exp().matrix().asDiagonal() * C.adjoint();

        // store the propagator matrices
        R.at(i) = RI, K.at(i) = KI.asDiagonal();
    }

    // return the propagator
    return {R, K};
}

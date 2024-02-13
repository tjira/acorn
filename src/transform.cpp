#include "transform.h"

Matrix<> Transform::Single(const Matrix<>& Aao, const Matrix<>& C) {
    // create the transformed matrix
    Matrix<> Amo(Aao.rows(), Aao.cols());

    // perform the transformation
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < Aao.rows(); i++) {
        for (int j = 0; j < Aao.cols(); j++) {
            for (int a = 0; a < Amo.rows(); a++) {
                for (int b = 0; b < Amo.cols(); b++) {
                    Amo(i, j) += Aao(a, b) * C(a, i) * C(b, j);
                }
            }
        }
    }

    // return the transformed matrix
    return Amo;
}

Matrix<> Transform::SingleSpin(const Matrix<>& Aao, const Matrix<>& C) {
    // expand the dimentsions of the MO basis matrix
    Matrix<> Ams = Numpy::Repeat(Numpy::Repeat(Transform::Single(Aao, C), 2, 0), 2, 1);

    // create repeating zeros and ones and tile it to a matrix
    Vector<double> ind(Ams.rows()); std::iota(ind.begin(), ind.end(), 0); ind = ind.unaryExpr([](auto s) {return double(int(s) % 2);});
    Matrix<bool> indm(ind.size(), ind.size()); for (long i = 0; i < ind.size(); i++) indm.row(i) = ind.transpose().array() == ind(i);

    // element wise multiply and return
    return Ams.array() * indm.cast<double>().array();
}

Tensor<> Transform::Coulomb(const Tensor<>& Jao, const Matrix<>& Ca, const Matrix<>& Cb) {
    // declare the Coulomb integral in molecular orbital basis and tensors of partial transform
    Tensor<4> J01(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); J01.setZero();
    Tensor<4> J02(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); J02.setZero();
    Tensor<4> J03(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); J03.setZero();
    Tensor<4> Jmo(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); Jmo.setZero();

    // perform the transformation
    J01 = toTensor(Ca).contract(Jao, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 0)});
    J02 = toTensor(Ca).contract(J01, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 1)});
    J03 = toTensor(Cb).contract(J02, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 2)});
    Jmo = toTensor(Cb).contract(J03, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 3)});

    // return JMO
    return Jmo;
}

Tensor<> Transform::Coulomb(const Tensor<>& Jao, const Matrix<>& C) {
    return Coulomb(Jao, C, C);
}

Tensor<> Transform::CoulombSpin(const Tensor<>& Jao, const Matrix<>& C) {
    // create the spin indices
    Vector<> ind(2 * C.rows()); std::iota(ind.begin(), ind.end(), 0); ind = ind.unaryExpr([](auto s) {return double(int(s) % 2);});
    Matrix<bool> indm(ind.size(), ind.size()); for (int i = 0; i < ind.size(); i++) indm.row(i) = ind.transpose().array() == ind(i);

    // expansion coefficients in MS basis and the Coulomb tensor in AS basis
    Matrix<> Cs = Numpy::Repeat(C, 2, 1).replicate<2, 1>().array() * Numpy::Repeat(indm.cast<double>(), C.rows(), 0).topRows(indm.rows()).array();
    Tensor<4> Js = Numpy::Kron(Matrix<>::Identity(2, 2), Numpy::Kron(Matrix<>::Identity(2, 2), Jao).shuffle(Eigen::array<int, 4>{3, 2, 1, 0}));

    // return coulomb tensor in MS
    return Coulomb(Js, Cs);
}

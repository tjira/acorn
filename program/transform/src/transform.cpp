#include "transform.h"

Tensor<> Transform::CoulombSpatial(const Tensor<>& Jao, const Matrix<>& Cmo) {
    // perform the transform
    Tensor<4> J01 = Cmo.tensor().contract<4, 1>(Jao, {Eigen::IndexPair<int>(0, 0)});
    Tensor<4> J02 = Cmo.tensor().contract<4, 1>(J01, {Eigen::IndexPair<int>(0, 1)});
    Tensor<4> J03 = Cmo.tensor().contract<4, 1>(J02, {Eigen::IndexPair<int>(0, 2)});
    Tensor<4> Jmo = Cmo.tensor().contract<4, 1>(J03, {Eigen::IndexPair<int>(0, 3)});

    return Jmo; // return the Coulomb integrals in molecular orbital basis
}

Tensor<> Transform::CoulombSpin(const Tensor<>& Jao, const Matrix<>& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times
    Matrix<> P = Matrix<>::Identity(Cmo.rows()).repeat(2, 1);

    // create the spin masks
    Matrix<> M(Cmo.rows(), 2 * Cmo.cols(), [](int, int j) {return 1 - j % 2;});
    Matrix<> N(Cmo.rows(), 2 * Cmo.cols(), [](int, int j) {return j % 2;});

    // transform the wfn coefficients to the spin basis
    Matrix<> Cms = Cmo.dot(P).vjoin(Cmo.dot(P)) * M.vjoin(N);

    // return the transformed matrix
    return CoulombSpatial(Matrix<>::Identity(2).kron(Matrix<>::Identity(2).kron(Jao).t({3, 2, 1, 0})), Cms);
}

Matrix<> Transform::SingleSpatial(const Matrix<>& Aao, const Matrix<>& Cmo) {
    // perform the transform
    Matrix<> A01 = Cmo.tensor().contract<2, 1>(Aao.tensor(), {Eigen::IndexPair<int>(0, 0)}).matrix();
    Matrix<> Amo = Cmo.tensor().contract<2, 1>(A01.tensor(), {Eigen::IndexPair<int>(0, 1)}).matrix();

    return Amo; // return the single-electron integrals in molecular orbital basis
}

Matrix<> Transform::SingleSpin(const Matrix<>& Aao, const Matrix<>& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times
    Matrix<> P = Matrix<>::Identity(Cmo.rows()).repeat(2, 1);

    // create the spin masks
    Matrix<> M(Cmo.rows(), 2 * Cmo.cols(), [](int, int j) {return 1 - j % 2;});
    Matrix<> N(Cmo.rows(), 2 * Cmo.cols(), [](int, int j) {return j % 2;});

    // transform the wfn coefficients to the spin basis
    Matrix<> Cms = Cmo.dot(P).vjoin(Cmo.dot(P)) * M.vjoin(N);

    // return the transformed matrix
    return SingleSpatial(Matrix<>::Identity(2).kron(Aao), Cms);
}

#include "transform.h"

Tensor<4> Transform::CoulombSpatial(const Tensor<4>& Jao, Matrix& Cmo) {
    // perform the transform
    Tensor<4> J01 = TENSORMAP(Cmo).contract(Jao, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 0)});
    Tensor<4> J02 = TENSORMAP(Cmo).contract(J01, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 1)});
    Tensor<4> J03 = TENSORMAP(Cmo).contract(J02, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 2)});
    Tensor<4> Jmo = TENSORMAP(Cmo).contract(J03, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 3)});

    // return the Coulomb integrals in molecular orbital basis
    return Jmo;
}

Tensor<4> Transform::CoulombSpin(const Tensor<4>& Jao, const Matrix& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times and define the spin mask
    Matrix P = Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.rows(), std::function<double(int, int)>([](int i, int j) {return j == 2 * i || j == 2 * i + 1;})), MN(2 * Cmo.cols(), 2 * Cmo.cols());

    // initialize the spin mask
    MN << Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 1 - j % 2;})),
          Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 0 + j % 2;}));

    // transform the wfn coefficients to the spin basis
    Matrix Cms = (Cmo * P).replicate<2, 1>().cwiseProduct(MN);

    // return the transformed matrix
    return CoulombSpatial(Eigen::Kron(Matrix::Identity(2, 2), Eigen::Kron(Matrix::Identity(2, 2), Jao).shuffle(Eigen::array<int, 4>{3, 2, 1, 0})), Cms);
}

Matrix Transform::SingleSpatial(const Matrix& Aao, Matrix& Cmo) {return Cmo.transpose() * Aao * Cmo;}

Matrix Transform::SingleSpin(const Matrix& Aao, const Matrix& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times and define the spin mask
    Matrix P = Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.rows(), std::function<double(int, int)>([](int i, int j) {return j == 2 * i || j == 2 * i + 1;})), MN(2 * Cmo.cols(), 2 * Cmo.cols());

    // initialize the spin mask
    MN << Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 1 - j % 2;})),
          Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 0 + j % 2;}));

    // transform the wfn coefficients to the spin basis
    Matrix Cms = (Cmo * P).replicate<2, 1>().cwiseProduct(MN);

    // return the transformed matrix
    return SingleSpatial(Eigen::kroneckerProduct(Matrix::Identity(2, 2), Aao), Cms);
}

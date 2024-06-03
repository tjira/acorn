#include "transform.h"

EigenTensor<> Transform::CoulombSpatial(const EigenTensor<>& Jao, EigenMatrix<>& Cmo) {
    // perform the transform
    EigenTensor<> J01 = TENSORMAP(Cmo).contract(Jao, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 0)});
    EigenTensor<> J02 = TENSORMAP(Cmo).contract(J01, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 1)});
    EigenTensor<> J03 = TENSORMAP(Cmo).contract(J02, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 2)});
    EigenTensor<> Jmo = TENSORMAP(Cmo).contract(J03, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 3)});

    // return the Coulomb integrals in molecular orbital basis
    return Jmo;
}

EigenTensor<> Transform::CoulombSpin(const EigenTensor<>& Jao, const EigenMatrix<>& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times and define the spin mask
    EigenMatrix<> P = EigenMatrix<>::NullaryExpr(Cmo.rows(), 2 * Cmo.rows(), std::function<double(int, int)>([](int i, int j) {return j == 2 * i || j == 2 * i + 1;})), MN(2 * Cmo.cols(), 2 * Cmo.cols());

    // initialize the spin mask
    MN << EigenMatrix<>::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 1 - j % 2;})),
          EigenMatrix<>::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 0 + j % 2;}));

    // transform the wfn coefficients to the spin basis
    EigenMatrix<> Cms = (Cmo * P).replicate<2, 1>().cwiseProduct(MN);

    // return the transformed matrix
    return CoulombSpatial(Eigen::Kron<double>(EigenMatrix<>::Identity(2, 2), Eigen::Kron<double>(EigenMatrix<>::Identity(2, 2), Jao).shuffle(Eigen::array<int, 4>{3, 2, 1, 0})), Cms);
}

EigenMatrix<> Transform::SingleSpatial(const EigenMatrix<>& Aao, EigenMatrix<>& Cmo) {return Cmo.transpose() * Aao * Cmo;}

EigenMatrix<> Transform::SingleSpin(const EigenMatrix<>& Aao, const EigenMatrix<>& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times and define the spin mask
    EigenMatrix<> P = EigenMatrix<>::NullaryExpr(Cmo.rows(), 2 * Cmo.rows(), std::function<double(int, int)>([](int i, int j) {return j == 2 * i || j == 2 * i + 1;})), MN(2 * Cmo.cols(), 2 * Cmo.cols());

    // initialize the spin mask
    MN << EigenMatrix<>::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 1 - j % 2;})),
          EigenMatrix<>::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 0 + j % 2;}));

    // transform the wfn coefficients to the spin basis
    EigenMatrix<> Cms = (Cmo * P).replicate<2, 1>().cwiseProduct(MN);

    // return the transformed matrix
    return SingleSpatial(Eigen::kroneckerProduct(EigenMatrix<>::Identity(2, 2), Aao), Cms);
}

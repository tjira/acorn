#include "transform.h"

Tensor<> Transform::CoulombSpatial(const Tensor<>& Jao, const Matrix<>& Cmo) {
    // declare the Coulomb integral in molecular orbital basis and tensors of partial transform
    Tensor<4> J01(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); J01.zero();
    Tensor<4> J02(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); J02.zero();
    Tensor<4> J03(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); J03.zero();
    Tensor<4> Jmo(Jao.dimension(0), Jao.dimension(1), Jao.dimension(2), Jao.dimension(3)); Jmo.zero();

    // perform the transformation
    J01 = Cmo.tensor().contract(Jao, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 0)});
    J02 = Cmo.tensor().contract(J01, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 1)});
    J03 = Cmo.tensor().contract(J02, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 2)});
    Jmo = Cmo.tensor().contract(J03, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 3)});

    return Jmo; // return the Coulomb integral in molecular orbital basis
}

Tensor<> Transform::CoulombSpin(const Tensor<>& Jao, const Matrix<>& Cmo) {
    // create the spin indices
    // Matrix<> ind(2 * Cmo.rows(), 1); std::iota(ind.begin(), ind.end(), 0);// ind = ind.unaryExpr([](auto s) {return double(int(s) % 2);});
    // Matrix<bool> indm(ind.size(), ind.size()); for (int i = 0; i < ind.size(); i++) indm.row(i) = ind.transpose().array() == ind(i);

    // expansion coefficients in MS basis and the Coulomb tensor in AS basis
    // Matrix<> Cs = Numpy::Repeat(C, 2, 1).replicate<2, 1>().array() * Numpy::Repeat(indm.cast<double>(), C.rows(), 0).topRows(indm.rows()).array();
    // Tensor<4> Js = Numpy::Kron(Matrix<>::Identity(2, 2), Numpy::Kron(Matrix<>::Identity(2, 2), Jao).shuffle(Eigen::array<int, 4>{3, 2, 1, 0}));

    // return coulomb tensor in MS
    // return Coulomb(Js, Cs);
}

#include "transform.h"

Tensor<> Transform::Coulomb(const Tensor<>& Jao, const Matrix<>& Cmo) {
    // perform the 4 partial transforms of 5th order
    Tensor<> J01 = Cmo.tensor().contract(Jao, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 0)});
    Tensor<> J02 = Cmo.tensor().contract(J01, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 1)});
    Tensor<> J03 = Cmo.tensor().contract(J02, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 2)});
    Tensor<> Jmo = Cmo.tensor().contract(J03, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 3)});

    return Jmo; // return the Coulomb integral in molecular orbital basis
}

#pragma once

#include "numpy.h"

namespace Transform {
    // single electron transforms
    Matrix<> SingleSpin(const Matrix<>& A, const Matrix<>& C);
    Matrix<> Single(const Matrix<>& A, const Matrix<>& C);

    // double electron transforms
    Tensor<> Coulomb(const Tensor<>& J, const Matrix<>& C);
}

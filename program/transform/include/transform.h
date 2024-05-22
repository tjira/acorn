#pragma once

#include "tensor.h"

namespace Transform {
    // two electron integrals
    Tensor<> CoulombSpatial(const Tensor<>& Jao, const Matrix<>& Cmo);
    Tensor<> CoulombSpin(const Tensor<>& Jao, const Matrix<>& Cmo);

    // one electron integrals
    Matrix<> SingleSpatial(const Matrix<>& Aao, const Matrix<>& Cmo);
    Matrix<> SingleSpin(const Matrix<>& Aao, const Matrix<>& Cmo);
}

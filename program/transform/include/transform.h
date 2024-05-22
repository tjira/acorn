#pragma once

#include "tensor.h"

namespace Transform {
    Tensor<> CoulombSpatial(const Tensor<>& Jao, const Matrix<>& Cmo);
    Tensor<> CoulombSpin(const Tensor<>& Jao, const Matrix<>& Cmo);
}

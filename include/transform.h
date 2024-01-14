#pragma once

#include "default.h"
#include "eigen.h"

namespace Transform {
    Tensor<> Coulomb(const Tensor<>& J, const Matrix<>& C);
}

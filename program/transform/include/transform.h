#pragma once

#include "tensor.h"

namespace Transform {
    Tensor<> Coulomb(const Tensor<>& Jao, const Matrix<>& Cmo);
}

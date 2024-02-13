#pragma once

#include "eigen.h"

namespace Numpy {
    std::vector<std::vector<int>> Combinations(int n, int k);
    Matrix<> Repeat(const Matrix<>& A, int count, int axis);
    Tensor<> Kron(const Matrix<>& A, const Tensor<>& B);
}

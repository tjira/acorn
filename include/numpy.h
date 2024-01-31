#pragma once

#include "eigen.h"

namespace Numpy {
    // useful math functions
    std::vector<std::vector<int>> Combinations(int n, int k);
    Matrix<> Repeat(const Matrix<>& A, int count, int axis);
}

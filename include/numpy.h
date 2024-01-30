#pragma once

#include "eigen.h"

namespace Numpy {
    // useful math functions
    std::vector<std::vector<int>> Combinations(int n, int k);
    Matrix<> Repeat(const Matrix<>& A, int count, int axis);
    long Factorial(long n);

    // templated functions
    template <size_t D> Tensor<D> Kron(const Matrix<>& A, const Tensor<D>& B);
}

template <size_t D>
Tensor<D> Numpy::Kron(const Matrix<>& A, const Tensor<D>& B) {
    // define the resulting tensor
    Tensor<D> C(B.dimension(0), B.dimension(1), A.rows() * B.dimension(2), A.cols() * B.dimension(3));

    // perform the cronecker product
    for (int i = 0; i < C.dimension(0); i++) {
        for (int j = 0; j < C.dimension(1); j++) {
            for (int k = 0; k < C.dimension(2); k++) {
                for (int l = 0; l < C.dimension(3); l++) {
                    C(i, j, k, l) = A(k / B.dimension(2), l / B.dimension(3)) * B(i, j, k % B.dimension(2), l % B.dimension(3));
                }
            }
        }
    }

    // return the result
    return C;
}

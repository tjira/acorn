#pragma once

#include "eigen.h"

namespace Numpy {
    Vector<std::complex<double>> FFT(Vector<std::complex<double>> in, int sign);
    Vector<std::complex<double>> IFFT(Vector<std::complex<double>> in);
    Vector<std::complex<double>> FFT(Vector<std::complex<double>> in);
    std::vector<std::vector<int>> Combinations(int n, int k);
    Matrix<> Repeat(const Matrix<>& A, int count, int axis);
    Tensor<> Kron(const Matrix<>& A, const Tensor<>& B);
}

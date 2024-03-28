#pragma once

#include "eigen.h"

namespace Numpy {
    void FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int> sizes, int sign);
    Matrix<std::complex<double>> FFT(Matrix<std::complex<double>> in, int sign = -1);
    double Moment(const Vector<>& x, const Vector<>& f, double c, int n);
    void HFFT1D(std::complex<double>* in, double* out, int size);
    std::vector<std::vector<int>> Combinations(int n, int k);
    Matrix<> Repeat(const Matrix<>& A, int count, int axis);
    Vector<double> HFFT1D(Vector<std::complex<double>> in);
    Tensor<> Kron(const Matrix<>& A, const Tensor<>& B);
}

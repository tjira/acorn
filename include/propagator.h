#pragma once

#include "matrixofmatrices.h"
#include "modelsystem.h"

template <int S>
class Propagator {
public:
    static Matrix<std::complex<double>> Propagate(const MatrixOfMatrices<Matrix<std::complex<double>>>& R, const MatrixOfMatrices<Matrix<std::complex<double>>>& K, Matrix<std::complex<double>> psi);
    static std::tuple<ComplexMatrixOfMatrices, ComplexMatrixOfMatrices> Get(const ModelSystem& system, const ComplexMatrixOfMatrices& V, const Matrix<>& ksq, double step, bool real);
    static double Energy(const ModelSystem& system, const MatrixOfMatrices<Matrix<>>& V, const Matrix<>& r, const Matrix<>& ksq, Matrix<std::complex<double>> psi);
};

#pragma once

#define MATRIXMAP(A) Eigen::      Map<Matrix   >(A.data(), A.dimension(0), A.dimension(1))
#define TENSORMAP(A) Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(),       A.cols()      )

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

// general templated eigen types
template <typename T, size_t D> using EigenTensor = Eigen::Tensor<T, D,                                 Eigen::ColMajor>;
template <typename T>           using EigenMatrix = Eigen::Matrix<T,    Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <typename T>           using EigenVector = Eigen::Matrix<T,    Eigen::Dynamic, 1             , Eigen::ColMajor>;


// specific eigen types
template <size_t D> using Tensor = EigenTensor<double, D>; template <size_t D> using IntegerTensor = EigenTensor<int, D>; template <size_t D> using ComplexTensor = EigenTensor<std::complex<double>, D>;
                    using Matrix = EigenMatrix<double   >;                     using IntegerMatrix = EigenMatrix<int   >;                     using ComplexMatrix = EigenMatrix<std::complex<double>   >;
                    using Vector = EigenVector<double   >;                     using IntegerVector = EigenVector<int   >;                     using ComplexVector = EigenVector<std::complex<double>   >;

namespace Eigen {
    // general functions
    EigenTensor<double, 4> Kron(const EigenMatrix<double>& A, const EigenTensor<double, 4>& B);

    // file readers
    EigenTensor<double, 4> LoadTensor(const std::string& path);
    EigenMatrix<double   > LoadMatrix(const std::string& path);

    // file writers
    void Write(const std::string& path, const EigenTensor<double, 4>& A);
    void Write(const std::string& path, const EigenMatrix<double   >& A);
}

// operators
std::ostream& operator<<(std::ostream& os, const Matrix& A);

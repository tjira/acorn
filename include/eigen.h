#pragma once

#include "constant.h"

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

// include the fourier transform
#include <unsupported/Eigen/FFT>

// define the Eigen types
template <typename T = double> using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <size_t D = 4, typename T = double> using Tensor = Eigen::Tensor<T, D, Eigen::ColMajor>;
template <typename T = double> using Vector = Eigen::Vector<T, Eigen::Dynamic>;

// conversion functions
inline Matrix<> toMatrix(Tensor<2> A) {return Eigen::Map<Matrix<>>(A.data(), A.dimension(0), A.dimension(1));}
inline Tensor<2> toTensor(Matrix<> A) {return Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(), A.cols());}

// printing functions
template <typename T> std::ostream& operator<<(std::ostream& os, const Tensor<4, T>& A);
template <typename T> std::ostream& operator<<(std::ostream& os, const Matrix<T>& A);
template <typename T> std::ostream& operator<<(std::ostream& os, const Vector<T>& A);

// exporting functions
template <typename T> void EigenWrite(const std::string& fname, const Tensor<4, T>& A);
template <typename T> void EigenWrite(const std::string& fname, const Matrix<T>& A);
template <typename T> void EigenWrite(const std::string& fname, const Vector<T>& A);

//converter functions
template <typename T> T EigenConj(const T& A) {return A.unaryExpr([](auto x) {return std::conj(x);});}

// function to fourier transform a complex vector
inline Vector<std::complex<double>> EigenFourier(const Vector<std::complex<double>>& A, bool inv = false) {Eigen::FFT<double> fft; Vector<std::complex<double>> B(A.size()); inv ? fft.inv(B, A) : fft.fwd(B, A); return B;}

#include <fstream>
#include <iomanip>

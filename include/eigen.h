#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

template <typename T = double> using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <size_t D = 4, typename T = double> using Tensor = Eigen::Tensor<T, D, Eigen::ColMajor>;
template <typename T = double> using Vector = Eigen::Vector<T, Eigen::Dynamic>;

inline Matrix<> toMatrix(Tensor<2> A) {return Eigen::Map<Matrix<>>(A.data(), A.dimension(0), A.dimension(1));}
inline Tensor<2> toTensor(Matrix<> A) {return Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(), A.cols());}

#pragma once

#include "constant.h"

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

// eigen include
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop

// include necessities
#include <filesystem>

// define the Eigen types
template <typename T = double> using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <size_t D = 4, typename T = double> using Tensor = Eigen::Tensor<T, D, Eigen::ColMajor>;
template <typename T = double> using Vector = Eigen::Vector<T, Eigen::Dynamic>;

// conversion functions
inline Matrix<> toMatrix(Tensor<2> A) {return Eigen::Map<Matrix<>>(A.data(), A.dimension(0), A.dimension(1));}
inline Tensor<2> toTensor(Matrix<> A) {return Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(), A.cols());}
inline Vector<> toVector(Tensor<1> A) {return Eigen::Map<Vector<>>(A.data(), A.dimension(0));}

// printing functions
template <typename T> std::ostream& operator<<(std::ostream& os, const Tensor<3, T>& A);
template <typename T> std::ostream& operator<<(std::ostream& os, const Tensor<4, T>& A);
template <typename T> std::ostream& operator<<(std::ostream& os, const Tensor<5, T>& A);
template <typename T> std::ostream& operator<<(std::ostream& os, const Matrix<T>& A);
template <typename T> std::ostream& operator<<(std::ostream& os, const Vector<T>& A);

// exporting functions
template <typename T> void EigenWrite(const std::filesystem::path& fname, const Tensor<3, T>& A);
template <typename T> void EigenWrite(const std::filesystem::path& fname, const Tensor<4, T>& A);
template <typename T> void EigenWrite(const std::filesystem::path& fname, const Tensor<5, T>& A);
template <typename T> void EigenWrite(const std::filesystem::path& fname, const Matrix<T>& A);
template <typename T> void EigenWrite(const std::filesystem::path& fname, const Vector<T>& A);

//converter functions
template <typename T> T EigenConj(const T& A) {return A.unaryExpr([](auto x) {return std::conj(x);});}

// complamentary functions
inline bool VectorContains(const std::vector<int>& v, const int& e) {return std::find(v.begin(), v.end(), e) != v.end();}
inline bool StringContains(const std::string& s, const char& e) {return s.find(e) != std::string::npos;}

// number of threads and input path global variable
inline std::filesystem::path ip; inline int nthread = 1;

// include necessities
#include <fstream>
#include <iomanip>

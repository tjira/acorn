#pragma once

#define MATRIXMAP(A) Eigen::Map<EigenMatrix<>>(A.data(), A.dimension(0), A.dimension(1))
#define TENSORMAP(A) Eigen::TensorMap<EigenTensor<2>>(A.data(), A.rows(), A.cols())

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

#include <fstream>

template <typename T = double> using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <size_t D = 4, typename T = double> using EigenTensor = Eigen::Tensor<T, D, Eigen::ColMajor>;
template <typename T = double> using EigenVector = Eigen::Vector<T, Eigen::Dynamic>;

namespace Eigen {
    // complementary constructors
    template <typename T = double> EigenMatrix<T> IndexFunction(int rows, int cols, const std::function<T(int, int)>& func);

    // general functions
    template <typename T = double> EigenTensor<4, T> Kron(const EigenMatrix<T>& A, const EigenTensor<4, T>& B);
    template <typename T = double> EigenMatrix<   T> Kron(const EigenMatrix<T>& A, const EigenMatrix<   T>& B);
    template <typename T = double> EigenMatrix<T> Repeat(const EigenMatrix<T>& A, int count, int axis);

    // file readers
    template <typename T = double> EigenTensor<4, T> LoadTensor(const std::string& path);
    template <typename T = double> EigenMatrix<   T> LoadMatrix(const std::string& path);

    // file writers
    template <typename T = double> void Write(const std::string& path, const EigenTensor<4, T>& A);
    template <typename T = double> void Write(const std::string& path, const EigenMatrix<   T>& A);
}

template <typename T>
EigenMatrix<T> Eigen::IndexFunction(int m, int n, const std::function<T(int, int)>& func) {
    EigenMatrix<T> A(m, n); for (int i = 0; i < m; i++) {for (int j = 0; j < n; j++) A(i, j) = func(i, j);} return A;
}

template <typename T>
EigenMatrix<T> Eigen::Kron(const EigenMatrix<T>& A, const EigenMatrix<T>& B) {
    // define the matrix where the product will be stored
    EigenMatrix<T> C(A.rows() * B.rows(), A.cols() * B.cols());

    // perform the Kronecker product
    for (int i = 0; i < C.rows(); i++) {
        for (int j = 0; j < C.cols(); j++) {
            C(i, j) = A(i / B.rows(), j / B.cols()) * B(i % B.rows(), j % B.cols());
        }
    }

    // return the Kronecker product
    return C;
}

template <typename T>
EigenTensor<4, T> Eigen::Kron(const EigenMatrix<T>& A, const EigenTensor<4, T>& B) {
    // define the tensor where the product will be stored
    EigenTensor<4, T> C(B.dimension(0), B.dimension(1), A.rows() * B.dimension(2), A.cols() * B.dimension(3));

    // perform the Kronecker product
    for (int i = 0; i < C.dimension(0); i++) {
        for (int j = 0; j < C.dimension(1); j++) {
            for (int k = 0; k < C.dimension(2); k++) {
                for (int l = 0; l < C.dimension(3); l++) {
                    C(i, j, k, l) = A(k / B.dimension(2), l / B.dimension(3)) * B(i, j, k % B.dimension(2), l % B.dimension(3));
                }
            }
        }
    }

    // return the Kronecker product
    return C;
}

template <typename T>
EigenMatrix<T> Eigen::Repeat(const EigenMatrix<T>& A, int count, int axis) {
    // create the new matrix with the repeated dimensions
    EigenMatrix<> B(axis ? A.rows() : count * A.rows(), axis ? count * A.cols() : A.cols());

    // throw error if axis is out of bounds
    if (axis != 0 && axis != 1) throw std::runtime_error("UNKNOWN AXIS IN MATRIX REPEAT");

    // repeat rows for axis 0 and columns for axis 1
    if (axis == 0) for (int i = 0; i < A.rows(); i++) for (int j = 0; j < count; j++) B.row(i * count + j) = A.row(i);
    else if (axis == 1) for (int i = 0; i < A.cols(); i++) for (int j = 0; j < count; j++) B.col(i * count + j) = A.col(i);

    // return the repeated matrix
    return B;
}

template <typename T>
EigenMatrix<T> Eigen::LoadMatrix(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and create the tensor
    int rows, cols; file >> rows >> cols; EigenMatrix<T> A(rows, cols);

    // read the tensor by dimensions, assign the values and return the matrix
    for (int i = 0; i < rows; i++) {for (int j = 0; j < cols; j++) file >> A(i, j);} return A;
}

template <typename T>
EigenTensor<4, T> Eigen::LoadTensor(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and assign them to an array
    std::array<int, 4> dims; for (int i = 0; i < 4; i++) file >> dims.at(i);

    // create the tensor
    EigenTensor<4, T> A(dims.at(0), dims.at(1), dims.at(2), dims.at(3));

    // read the tensor by dimensions, assign the values and return the tensor
    for (int i = 0; i < dims.at(0); i++) for (int j = 0; j < dims.at(1); j++) {
        for (int k = 0; k < dims.at(2); k++) for (int l = 0; l < dims.at(3); l++) file >> A(i, j, k, l);
    }

    // return the tensor
    return A;
}

template <typename T>
void Eigen::Write(const std::string& path, const EigenTensor<4, T>& A) {
    // open the output file and extract the dimensions
    std::ofstream file(path); auto dims = A.dimensions();

    // write the dimensions to the header and set the formatting
    for (size_t i = 0; i < dims.size(); i++) {file << dims.at(i) << (i < dims.size() - 1 ? " " : "");} file << "\n" << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < dims.at(0); i++) {
        for (int j = 0; j < dims.at(1); j++) {
            for (int k = 0; k < dims.at(2); k++) {
                for (int l = 0; l < dims.at(3); l++) {
                    file << std::setw(20) << A(i, j, k, l) << (k < dims.at(2) - 1 ? " " : "");
                }
            } file << "\n";
        }
    }
}

template <typename T>
void Eigen::Write(const std::string& path, const EigenMatrix<T>& A) {
    // open the output file and write the dimensions to the header
    std::ofstream file(path); file << A.rows() << " " << A.cols() << "\n" << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < A.rows(); i++, file << "\n") for (int j = 0; j < A.cols(); j++) file << std::setw(20) << A(i, j) << (j < A.cols() - 1 ? " " : "");
}

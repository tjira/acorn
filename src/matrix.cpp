#include "matrix.h"

#include <fstream>

// static functions
template <typename T> Matrix<T> Matrix<T>::Load(const std::string& path) {return Tensor<2, T>::Load(path).matrix();}
template <typename T> Matrix<T> Matrix<T>::Identity(int n) {return Matrix<T>(EigenMatrix<T>::Identity(n, n));}

// matrix operators
template <typename T> Matrix<T> Matrix<T>::operator*(const Matrix<T>& A) const {return Matrix<T>(mat.cwiseProduct(A.mat));}
template <typename T> Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) const {return Matrix<T>(mat + A.mat);}
template <typename T> Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) const {return Matrix<T>(mat - A.mat);}

// special matrix operations
template <typename T> Matrix<T> Matrix<T>::dot(const Matrix<T>& A) const {return Matrix<T>(mat * A.mat);}
template <typename T> Matrix<T> Matrix<T>::t() const {return Matrix<T>(mat.transpose());}

// assignment operators and non-const functions
template <typename T> T& Matrix<T>::operator()(int i, int j) {return mat(i, j);}
template <typename T> T* Matrix<T>::data() {return mat.data();}
template <typename T> void Matrix<T>::zero() {mat.setZero();}

// block operations
template <typename T> Matrix<T> Matrix<T>::leftcols(int n) const {return Matrix<T>(mat.leftCols(n));}
template <typename T> const T& Matrix<T>::operator()(int i, int j) const {return mat(i, j);}
template <typename T> Matrix<T> Matrix<T>::col(int n) const {return Matrix<T>(mat.col(n));}
template <typename T> Matrix<T> Matrix<T>::row(int n) const {return Matrix<T>(mat.row(n));}
template <typename T> const T* Matrix<T>::data() const {return mat.data();}

// functions with number outputs
template <typename T> int Matrix<T>::cols() const {return mat.cols();}
template <typename T> int Matrix<T>::rows() const {return mat.rows();}
template <typename T> T Matrix<T>::norm() const {return mat.norm();}
template <typename T> T Matrix<T>::sum() const {return mat.sum();}

// input/output related functions
template <typename T> Tensor<2, T> Matrix<T>::tensor() const {EigenMatrix<T> temp = mat; return Tensor<2, T>(Eigen::TensorMap<EigenTensor<2, T>>(temp.data(), rows(), cols()));}

// self-adjoint eigenvalue problem solver
template <typename T> std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::eigh(const Matrix<T>& A) const {
    // solve the eigenvalue problem
    Eigen::GeneralizedSelfAdjointEigenSolver<EigenMatrix<T>> solver(mat, A.mat);

    // return the eigenvalues and eigenvectors
    return std::make_tuple(Matrix<T>(solver.eigenvalues()), Matrix<T>(solver.eigenvectors()));
}

// function to save a matrix to a file
template <typename T> void Matrix<T>::save(const std::string& path, std::vector<int> dims) const {
    // open the output file and set the dimensions if not provided
    std::ofstream file(path); if (dims.empty()) dims = {rows(), cols()};

    // write the dimensions to the header and set the formatting
    for (size_t i = 0; i < dims.size(); i++) {file << dims.at(i) << (i < dims.size() - 1 ? " " : "");} file << "\n" << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < rows(); i++, file << "\n") for (int j = 0; j < cols(); j++) file << std::setw(20) << mat(i, j) << (j < cols() - 1 ? " " : "");
}

template class Matrix<double>;

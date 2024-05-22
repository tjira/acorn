#include "matrix.h"

#include <fstream>

// constructors
template <typename T> Matrix<T>::Matrix(int m, int n, T value) : mat(m, n) {mat.fill(value);}

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

// constructor apply a function to each element of the matrix based on the indices
template <typename T> Matrix<T>::Matrix(int m, int n, const std::function<T(int, int)>& func) : mat(m, n) {
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) mat(i, j) = func(i, j);
}

// function to perform the Kronecker product with a matrix
template <typename T> Matrix<> Matrix<T>::kron(const Matrix<>& A) const {
    // define the matrix where the product will be stored
    Matrix<> newmat(rows() * A.rows(), cols() * A.cols());

    // perform the Kronecker product
    for (int i = 0; i < newmat.rows(); i++) {
        for (int j = 0; j < newmat.cols(); j++) {
            newmat(i, j) = mat(i / A.rows(), j / A.cols()) * A(i % A.rows(), j % A.cols());
        }
    }

    return newmat; // return the Kronecker product
}

// function to perform the Kronecker product with a 4th order tensor
template <typename T> Tensor<4, T> Matrix<T>::kron(const Tensor<4, T>& A) const {
    // define the tensor where the product will be stored
    Tensor<4, T> newten(A.dimension(0), A.dimension(1), rows() * A.dimension(2), cols() * A.dimension(3));

    // perform the Kronecker product
    for (int i = 0; i < newten.dimension(0); i++) {
        for (int j = 0; j < newten.dimension(1); j++) {
            for (int k = 0; k < newten.dimension(2); k++) {
                for (int l = 0; l < newten.dimension(3); l++) {
                    newten(i, j, k, l) = mat(k / A.dimension(2), l / A.dimension(3)) * A(i, j, k % A.dimension(2), l % A.dimension(3));
                }
            }
        }
    }

    return newten; // return the Kronecker product
}

// function to repeat the matrix along an axis
template <typename T> Matrix<T> Matrix<T>::repeat(int count, int axis) const {
    // create the new matrix with the repeated dimensions
    Matrix<> newmat(axis ? rows() : count * rows(), axis ? count * cols() : cols());

    // throw error if axis is out of bounds
    if (axis != 0 && axis != 1) throw std::runtime_error("UNKNOWN AXIS IN MATRIX REPEAT");

    // repeat rows for axis 0 and columns for axis 1
    if (axis == 0) for (int i = 0; i < rows(); i++) for (int j = 0; j < count; j++) newmat.mat.row(i * count + j) = mat.row(i);
    else if (axis == 1) for (int i = 0; i < cols(); i++) for (int j = 0; j < count; j++) newmat.mat.col(i * count + j) = mat.col(i);

    return newmat; // return the repeated matrix
}

// function to vertically join two matrices
template <typename T> Matrix<T> Matrix<T>::vjoin(const Matrix<T>& A) const {
    // create the new matrix with the joined dimensions
    Matrix<> newmat(rows() + A.rows(), cols());

    // copy the matrices into the new matrix
    newmat.mat.topRows(rows()) = mat; newmat.mat.bottomRows(A.rows()) = A.mat;

    return newmat; // return the vertically joined matrices
}

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

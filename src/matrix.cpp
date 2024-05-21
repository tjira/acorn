#include "matrix.h"

// matrix operators
template <typename T> Matrix<T> Matrix<T>::operator*(const Matrix<T>& A) const {return mat.cwiseProduct(A.mat);}
template <typename T> Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) const {return mat + A.mat;}
template <typename T> Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) const {return mat - A.mat;}

// special matrix operations
template <typename T> Matrix<T> Matrix<T>::dot(const Matrix<T>& A) const {return mat * A.mat;}
template <typename T> Matrix<T> Matrix<T>::t() const {return mat.transpose();}

// assignment operators
template <typename T> T& Matrix<T>::operator()(int i, int j) {return mat(i, j);}
template <typename T> T* Matrix<T>::data() {return mat.data();}

// block operations
template <typename T> Matrix<T> Matrix<T>::leftcols(int n) const {return mat.leftCols(n);}

// functions with number outputs
template <typename T> int Matrix<T>::cols() const {return mat.cols();}
template <typename T> int Matrix<T>::rows() const {return mat.rows();}
template <typename T> T Matrix<T>::norm() const {return mat.norm();}
template <typename T> T Matrix<T>::sum() const {return mat.sum();}

// non-const functions that alter the matrix values
template <typename T> void Matrix<T>::zero() {mat.setZero();}

// self-adjoint eigenvalue problem solver
template <typename T> std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::eigh(const Matrix<T>& A) const {
    Eigen::GeneralizedSelfAdjointEigenSolver<EigenMatrix<T>> solver(mat, A.mat); return {solver.eigenvalues(), solver.eigenvectors()};
}

// static function to load a matrix from a file
template <typename T> Matrix<T> Matrix<T>::Load(const std::string& filename) {
    std::ifstream file(filename); int rows, cols; file >> cols >> rows; Matrix<T> mat(rows, cols);
    for (int i = 0; i < rows; i++) {for (int j = 0; j < cols; j++) file >> mat(i, j);} return mat;
}

// function to save a matrix to a file
template <typename T> void Matrix<T>::save(const std::string& filename) const {
    std::ofstream file(filename); file << mat.rows() << " " << mat.cols() << std::endl << std::fixed << std::setprecision(14);
    for (int i = 0; i < rows(); i++, file << "\n") for (int j = 0; j < cols(); j++) file << std::setw(20) << mat(i, j) << " ";
}

// function to convert a matrix to a tensor
template <typename T> Tensor<2, T> Matrix<T>::tensor() const {
    EigenMatrix<T> temp = mat; return Eigen::TensorMap<EigenTensor<2>>(temp.data(), rows(), cols());
}

template class Matrix<double>;

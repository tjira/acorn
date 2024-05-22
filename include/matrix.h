#pragma once

#include <unsupported/Eigen/MatrixFunctions>

// define the Eigen matrix type
template <typename T = double> using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

// forward declare the tensor class
template <int D, typename T> class Tensor;

// define the matrix class
template <typename T = double> class Matrix {
    template <typename U> friend class Matrix;
public:
    Matrix(int m, int n) : mat(m, n) {mat.setZero();} Matrix(const EigenMatrix<T>& mat) : mat(mat) {}
    Matrix(int m, int n, const std::function<T(int, int)>& func); Matrix(int m, int n, T value);

    // static functions
    static Matrix<T> Load(const std::string& path); static Matrix<T> Identity(int n);

    // matrix operators
    Matrix<T> operator-(const Matrix<T>& A) const;
    Matrix<T> operator+(const Matrix<T>& A) const;
    Matrix<T> operator*(const Matrix<T>& A) const;

    // friend operators with scalars
    template <typename U> friend Matrix<U> operator*(double a, const Matrix<U>& A);
    template <typename U> friend Matrix<U> operator*(const Matrix<U>& A, double a);

    // special matrix operations
    Matrix<T> dot(const Matrix<T>& A) const; Matrix<T> t() const;
    Tensor<4, T> kron(const Tensor<4, T>& A) const;
    Matrix<> kron(const Matrix<>& A) const;

    // assignment operators and non-const functions
    T& operator()(int i, int j); T* data();

    // block operations
    const T& operator()(int i, int j) const; Matrix<T> leftcols(int n) const;
    Matrix<T> repeat(int count, int axis) const; const T* data() const;
    Matrix<> apply(const std::function<T(int, int)>& func) const;
    Matrix<T> row(int n) const; Matrix<T> col(int n) const;
    Matrix<T> vjoin(const Matrix<T>& A) const;

    // eigenproblem solvers
    std::tuple<Matrix<T>, Matrix<T>> eigh(const Matrix<T>& A) const;

    // functions with number outputs
    int cols() const; int rows() const; T norm() const; T sum() const;

    // input/output related functions
    void save(const std::string& path, std::vector<int> dims = {}) const; Tensor<2, T> tensor() const;
    template <typename U> friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& A);

private:
    EigenMatrix<T> mat;
};

// define the friend operators
template <typename T> std::ostream& operator<<(std::ostream& os, const Matrix<T>& A) {os << A.mat; return os;}
template <typename T> Matrix<T> operator*(double a, const Matrix<T>& A) {return Matrix<T>(a * A.mat);}
template <typename T> Matrix<T> operator*(const Matrix<T>& A, double a) {return Matrix<T>(a * A.mat);}

// define the tensor class
#include "tensor.h"

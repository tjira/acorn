#pragma once

#include <unsupported/Eigen/MatrixFunctions>

template <typename T = double> using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

template <int D, typename T> class Tensor; template <typename T = double> class Matrix {
    template <typename TT> friend class Matrix;
public:
    template <typename... dims> Matrix(dims... args);

    // static functions
    static Matrix<T> Load(const std::string& path);

    // matrix operators
    Matrix<T> operator-(const Matrix<T>& A) const;
    Matrix<T> operator+(const Matrix<T>& A) const;
    Matrix<T> operator*(const Matrix<T>& A) const;

    // friend operators with scalars
    template <typename TT> friend Matrix<TT> operator*(double a, const Matrix<TT>& A);
    template <typename TT> friend Matrix<TT> operator*(const Matrix<TT>& A, double a);

    // special matrix operations
    Matrix<T> dot(const Matrix<T>& A) const; Matrix<T> t() const;

    // assignment operators and non-const functions
    T& operator()(int i, int j); T* data(); void zero();

    // block operations
    Matrix<T> row(int n) const;  Matrix<T> col(int n) const; const T* data() const;
    const T& operator()(int i, int j) const; Matrix<T> leftcols(int n) const;

    // eigenproblem solvers
    std::tuple<Matrix<T>, Matrix<T>> eigh(const Matrix<T>& A) const;

    // functions with number outputs
    int cols() const; int rows() const; T norm() const; T sum() const;

    // input/output related functions
    void save(const std::string& path, std::vector<int> dims = {}) const; Tensor<2, T> tensor() const;

private:
    EigenMatrix<T> mat;
};

// the templated constructor with parameter pack
template <typename T> template <typename... aargs> Matrix<T>::Matrix(aargs... args) : mat(args...) {}

// define the friend operators
template <typename T> Matrix<T> operator*(double a, const Matrix<T>& A) {return a * A.mat;}
template <typename T> Matrix<T> operator*(const Matrix<T>& A, double a) {return a * A.mat;}

// define the tensor class
#include "tensor.h"

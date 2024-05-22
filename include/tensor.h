#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

// define the Eigen tensor type
template <size_t D = 4, typename T = double> using EigenTensor = Eigen::Tensor<T, D, Eigen::ColMajor>;

// define the integer concept
template <typename T> concept Integer = std::is_same_v<T, int>;

// forward declare the matrix class
template <typename T> class Matrix;

// define the tensor class
template <int D = 4, typename T = double> class Tensor {
    template <int E, typename U> friend class Tensor;
public:
    template <Integer... dims> Tensor(dims... args) : ten(args...) {} Tensor(const EigenTensor<D, T>& ten) : ten(ten) {}

    // static functions
    static Tensor<D, T> Load(const std::string& path);

    // tensor operators
    Tensor<D, T> operator+(const Tensor<D, T>& A) const;
    Tensor<D, T> operator-(const Tensor<D, T>& A) const;
    Tensor<D, T> operator*(const Tensor<D, T>& A) const;

    // friend operators with scalars
    template <int E, typename U> friend Tensor<E, U> operator*(double a, const Tensor<E, U>& A);
    template <int E, typename U> friend Tensor<E, U> operator*(const Tensor<E, U>& A, double a);

    // special tensor operations
    template <int E, long unsigned F> Tensor<D + E - 2 * F, T> contract(const Tensor<E, T>& A, const Eigen::array<Eigen::IndexPair<int>, F>& dims);
    Tensor<D, T> t(const Eigen::array<int, D>& axes) const;

    // assignment operators and non-const functions
    template <typename... dims> T& operator()(dims... args) {return ten(args...);} void zero();

    // block operations
    int dimension(int i) const;

    // input/output related functions
    void save(const std::string& path) const; Matrix<T> matrix() const;

private:
    EigenTensor<D, T> ten;
};

// define the friend operators
template <int D, typename T> Tensor<D, T> operator*(double a, const Tensor<D, T>& A) {return Tensor<D, T>(a * A.ten);}
template <int D, typename T> Tensor<D, T> operator*(const Tensor<D, T>& A, double a) {return Tensor<D, T>(a * A.ten);}

// function for tensor contraction
template <int D, typename T> template <int E, long unsigned F>
Tensor<D + E - 2 * F, T> Tensor<D, T>::contract(const Tensor<E, T>& A, const Eigen::array<Eigen::IndexPair<int>, F>& dims) {
    return Tensor<D + E - 2 * F, T>(ten.contract(A.ten, dims));
}

// define the matrix class
#include "matrix.h"

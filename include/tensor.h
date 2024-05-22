#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

template <size_t D = 4, typename T = double> using EigenTensor = Eigen::Tensor<T, D, Eigen::ColMajor>;

template <typename T> class Matrix; template <int D = 4, typename T = double> class Tensor {
    template <int DD, typename TT> friend class Tensor;
public:
    template <typename... aargs> Tensor(aargs... args);

    // static functions
    static Tensor<D, T> Load(const std::string& path);

    // tensor operators
    Tensor<D, T> operator+(const Tensor<D, T>& A) const;
    Tensor<D, T> operator-(const Tensor<D, T>& A) const;
    Tensor<D, T> operator*(const Tensor<D, T>& A) const;

    // friend operators with scalars
    template <int DD, typename TT> friend Tensor<DD, TT> operator*(double a, const Tensor<DD, TT>& A);
    template <int DD, typename TT> friend Tensor<DD, TT> operator*(const Tensor<DD, TT>& A, double a);

    // special matrix operations
    Tensor<2, T> contract(const Tensor<2, T>& A, const Eigen::array<Eigen::IndexPair<int>, 2>& dims);
    Tensor<D, T> t(const Eigen::array<int, D>& axes) const;

    // assignment operators and non-const functions
    template <typename... dims> T& operator()(dims... args) {return ten(args...);} void zero();

    // input/output related functions
    void save(const std::string& path) const; Matrix<T> matrix() const;

private:
    EigenTensor<D, T> ten;
};

// the templated constructor with parameter pack
template <int D, typename T> template <typename... aargs> Tensor<D, T>::Tensor(aargs... args) : ten(args...) {}

// define the friend operators
template <int D, typename T> Tensor<D, T> operator*(double a, const Tensor<D, T>& A) {return a * A.ten;}
template <int D, typename T> Tensor<D, T> operator*(const Tensor<D, T>& A, double a) {return a * A.ten;}

// define the matrix class
#include "matrix.h"

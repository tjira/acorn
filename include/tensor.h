#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

template <size_t D = 4, typename T = double> using EigenTensor = Eigen::Tensor<T, D, Eigen::ColMajor>;

template <typename T> class Matrix; template <int D = 4, typename T = double> class Tensor {
    template <int DD, typename TT> friend class Tensor;
public:
    template <typename... dims> Tensor(dims... args) : ten(args...) {}
    static Tensor<D, T> Load(const std::string& path);
    
    // operators
    template <typename... dims> T& operator()(dims... args) {return ten(args...);}

    template <int DD, typename TT> friend Tensor<DD, TT> operator*(double a, const Tensor<DD, TT>& A);
    template <int DD, typename TT> friend Tensor<DD, TT> operator*(const Tensor<DD, TT>& A, double a);
    Tensor<D, T> operator+(const Tensor<D, T>& A) const;
    Tensor<D, T> operator-(const Tensor<D, T>& A) const;
    Tensor<D, T> operator*(const Tensor<D, T>& A) const;

    // converters
    Matrix<T> matrix();
    Eigen::Tensor<T, D, Eigen::ColMajor> get() const {return ten;}
    void zero() {ten.setZero();}

    Tensor<2, T> contract(const Tensor<2, T>& A, const Eigen::array<Eigen::IndexPair<int>, 2>& dims);

    Tensor<D, T> transpose(const Eigen::array<int, D>& axes) const {return Tensor<D, T>(ten.shuffle(axes));}

    // input/output
    void save(const std::string& path) const;

private:
    EigenTensor<D, T> ten;
};

// define the friend operators
template <int D, typename T> Tensor<D, T> operator*(double a, const Tensor<D, T>& A) {return a * A.ten;}
template <int D, typename T> Tensor<D, T> operator*(const Tensor<D, T>& A, double a) {return a * A.ten;}

// define the matrix class
#include "matrix.h"

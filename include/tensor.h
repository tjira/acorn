#pragma once

#include "matrix.h"

template <int D = 4, typename T = double>
class Tensor {
public:
    template <typename... dims> Tensor(dims... args) : ten(args...) {}
    static Tensor<D, T> Load(const std::string& path);
    
    // operators
    template <typename... dims> T& operator()(dims... args) {return ten(args...);}

    Tensor<D, T> operator+(const Tensor<D, T>& other) const;
    Tensor<D, T> operator-(const Tensor<D, T>& other) const;
    Tensor<D, T> operator*(const T& a) const;

    // converters
    Matrix<T> matrix();
    Eigen::Tensor<T, D, Eigen::ColMajor> get() const {return ten;}
    void zero() {ten.setZero();}

    Tensor<2> contract(const Matrix<>& mat, const Eigen::array<Eigen::IndexPair<int>, 2>& dims);

    Tensor<D, T> transpose(const Eigen::array<int, D>& axes) const {return Tensor<D, T>(ten.shuffle(axes));}

    // input/output
    void save(const std::string& path) const;

    Eigen::Tensor<T, D, Eigen::ColMajor> ten;
};

#pragma once

#include <fstream>
#include <iomanip>

#include <unsupported/Eigen/CXX11/Tensor>

template <typename T = double>
class Matrix {
public:
    template <typename... dims> Matrix(dims... args) : mat(args...) {}
    static Matrix<T> Load(const std::string& path);

    // operators
    Matrix<T> operator-(const Matrix<T>& other) const;
    Matrix<T> operator+(const Matrix<T>& other) const;


    T& operator()(int i, int j) {return mat(i, j);}
    T& operator()(int i) {return mat.data()[i];}

    // getters
    int cols() const {return mat.cols();}
    int rows() const {return mat.rows();}
    void zero() {mat.setZero();}

    // output
    void save(const std::string& path) const;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> mat;
};

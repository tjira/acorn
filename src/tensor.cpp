#include "tensor.h"

#include <fstream>

// tensor operators
template <int D, typename T> Tensor<D, T> Tensor<D, T>::operator+(const Tensor<D, T>& A) const {return ten + A.ten;}
template <int D, typename T> Tensor<D, T> Tensor<D, T>::operator-(const Tensor<D, T>& A) const {return ten - A.ten;}
template <int D, typename T> Tensor<D, T> Tensor<D, T>::operator*(const Tensor<D, T>& A) const {return ten * A.ten;}

// special tensor operations
template <int D, typename T> Tensor<D, T> Tensor<D, T>::t(const Eigen::array<int, D>& axes) const {return ten.shuffle(axes);}

// block operations
template <int D, typename T> EigenTensor<D, T>::Dimensions Tensor<D, T>::dimensions() const {return ten.dimensions();}

// non-const functions
template <int D, typename T> void Tensor<D, T>::save(const std::string& path) const {return matrix().save(path, std::vector<int>(ten.dimensions().begin(), ten.dimensions().end()));}
template <int D, typename T> void Tensor<D, T>::zero() {ten.setZero();}

// function to convert a tensor to a matrix
template <int D, typename T>
Matrix<T> Tensor<D, T>::matrix() const {
    // get the dimensions of the resulting matrix
    int rows = 1, cols = 1; for (int i = 0; i < D; i++) if (i % 2) cols *= ten.dimension(i); else rows *= ten.dimension(i);

    // create the matrix and assign the values
    Matrix<T> mat(rows, cols); for (int i = 0; i < ten.size(); i++) {mat(i / cols, i % rows) = ten.data()[i];} return mat;
}

template <>
Tensor<2, double> Tensor<2, double>::Load(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and create the tensor
    int rows, cols; file >> cols >> rows; Tensor<2, double> ten(rows, cols);

    // read the tensor by dimensions, assign the values and return the tensor
    for (int i = 0; i < rows; i++) {for (int j = 0; j < cols; j++) file >> ten(i, j);} return ten;
}

template <>
Tensor<4, double> Tensor<4, double>::Load(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and assign them to an array
    std::array<int, 4> dims; for (int i = 0; i < 4; i++) file >> dims.at(i);

    // create the tensor
    Tensor<4, double> ten(dims.at(0), dims.at(1), dims.at(2), dims.at(3));

    // read the tensor by dimensions, assign the values and return the tensor
    for (int i = 0; i < dims.at(0); i++) for (int j = 0; j < dims.at(1); j++) {
        for (int k = 0; k < dims.at(2); k++) for (int l = 0; l < dims.at(3); l++) file >> ten(i, j, k, l);
    }

    return ten; // return the loaded tensor
}

template class Tensor<2, double>; template class Tensor<4, double>;

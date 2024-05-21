#include "tensor.h"

template <int D, typename T> Tensor<D, T> Tensor<D, T>::operator+(const Tensor<D, T>& A) const {return ten + A.ten;}
template <int D, typename T> Tensor<D, T> Tensor<D, T>::operator-(const Tensor<D, T>& A) const {return ten - A.ten;}
template <int D, typename T> Tensor<D, T> Tensor<D, T>::operator*(const Tensor<D, T>& A) const {return ten * A.ten;}

template <int D, typename T>
Tensor<2, T> Tensor<D, T>::contract(const Tensor<2, T>& A, const Eigen::array<Eigen::IndexPair<int>, 2>& dims) {
    return ten.contract(A.ten, dims);
}

template <int D, typename T>
Matrix<T> Tensor<D, T>::matrix() {
    if constexpr (D == 2) {
        Matrix<T> mat(ten.dimension(0), ten.dimension(1));
        for (int i = 0; i < ten.size(); i++) {mat(i / mat.cols(), i % mat.rows()) = ten.data()[i];} return mat;
    } else {


        Matrix<T> mat(ten.dimension(0) * ten.dimension(2), ten.dimension(1) * ten.dimension(3));
        for (int i = 0; i < ten.size(); i++) {mat(i / mat.cols(), i % mat.rows()) = ten.data()[i];} return mat;
    }
}

template <int D, typename T>
Tensor<D, T> Tensor<D, T>::Load(const std::string& path) {
    // open the input file and read the dimensions
    std::ifstream file(path); std::string line; std::getline(file, line);

    // read the tensor by dimensions
    if constexpr (D == 4) {
        Eigen::Tensor<T, D, Eigen::ColMajor> ten(7, 7, 7, 7);
        for (int i = 0; i < ten.size(); i++) file >> ten.data()[i];
        return Tensor<D, T>(ten);
    }
}

template <int D, typename T>
void Tensor<D, T>::save(const std::string& path) const {
    // open the output file
    std::ofstream file(path);

    if constexpr (D == 4) {
        // write the header and set the formatting
        file << ten.dimension(0) << " " << ten.dimension(1) << " " << ten.dimension(2) << " " << ten.dimension(3) << "\n" << std::fixed << std::setprecision(14);

        // write the matrix by rows
        for (int i = 0; i < ten.dimension(0); i++) {
            for (int j = 0; j < ten.dimension(1); j++) {
                for (int k = 0; k < ten.dimension(2); k++) {
                    for (int l = 0; l < ten.dimension(3); l++) {
                        file << std::setw(24) << ten(i, j, k, l) << " ";
                    }
                } file << "\n";
            }
        }
    }
}

template class Tensor<2, double>;
template class Tensor<4, double>;

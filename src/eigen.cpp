#include "eigen.h"

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& A) {
    // print the dimensions
    os << "(" << A.rows() << "x" << A.cols() << "):" << std::endl;

    // print the matrix
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            os << (j ? " " : "") << std::setw(20) << A(i, j);
        } os << (i < A.rows() - 1 ? "\n" : "");
    }

    // return stream
    return os;
}

template <typename T>
void EigenWrite(const std::string& fname, const Matrix<T>& A) {
    // open the output file and set the precision
    std::ofstream file(fname); file << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            file << std::setw(20) << A(i, j) << (j < A.cols() - 1 ? " " : "");
        } file << "\n";
    }
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Tensor<4, T>& A) {
    os << Matrix<T>(Eigen::Map<const Matrix<T>>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3))); return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& A) {
    os << Matrix<T>(Eigen::Map<const Matrix<T>>(A.data(), A.rows(), 1)); return os;
}

template <typename T>
void EigenWrite(const std::string& fname, const Tensor<4, T>& A) {
    EigenWrite(fname, Matrix<T>(Eigen::Map<const Matrix<T>>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3))));
}

template <typename T>
void EigenWrite(const std::string& fname, const Vector<T>& A) {
    EigenWrite(fname, Matrix<T>(Eigen::Map<const Matrix<T>>(A.data(), A.rows(), 1)));
}

template std::ostream& operator<<(std::ostream& os, const Tensor<4, double>& A);
template std::ostream& operator<<(std::ostream& os, const Matrix<double>& A);
template std::ostream& operator<<(std::ostream& os, const Vector<double>& A);

template void EigenWrite(const std::string& fname, const Tensor<4, double>& A);
template void EigenWrite(const std::string& fname, const Matrix<double>& A);
template void EigenWrite(const std::string& fname, const Vector<double>& A);

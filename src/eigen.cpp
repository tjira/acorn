#include "eigen.h"

template <typename T> std::ostream& operator<<(std::ostream& os, const Matrix<T>& A) {
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

template std::ostream& operator<<(std::ostream& os, const Matrix<double>& A);

#include "matrixofmatrices.h"

template <class M>
ComplexMatrixOfMatrices MatrixOfMatrices<M>::complex() const {
    // define the result
    ComplexMatrixOfMatrices result(data.size(), data.at(0).size());

    // fill the result with the results
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.at(i).size(); j++) {
            result(i, j) = data.at(i).at(j) * std::complex<double>(1, 0);
        }
    }

    // return the result
    return result;
}

template <class M>
RealMatrixOfMatrices MatrixOfMatrices<M>::real() const {
    // define the result
    RealMatrixOfMatrices result(data.size(), data.at(0).size());

    // fill the result with the results
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.at(i).size(); j++) {
            result(i, j) = data.at(i).at(j).real();
        }
    }

    // return the result
    return result;
}

template class MatrixOfMatrices<Eigen::Matrix<std::complex<double>, -1, -1, 0, -1, -1>>;
template class MatrixOfMatrices<Eigen::Matrix<double, -1, -1, 0, -1, -1>>;

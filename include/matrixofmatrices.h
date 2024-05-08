#pragma once

#include "eigen.h"
#include "numpy.h"

template <class M>
class MatrixOfMatrices {
public:
    MatrixOfMatrices(int rows, int cols) : data(rows, std::vector<M>(cols)) {}
    const M& operator()(int i, int j) const {return data.at(i).at(j);}
    MatrixOfMatrices<Matrix<std::complex<double>>> complex() const;
    M& operator()(int i, int j) {return data.at(i).at(j);}
    MatrixOfMatrices<Matrix<>> real() const;

private:
    std::vector<std::vector<M>> data;
};

typedef MatrixOfMatrices<Matrix<std::complex<double>>> ComplexMatrixOfMatrices;
typedef MatrixOfMatrices<Matrix<double>> RealMatrixOfMatrices;

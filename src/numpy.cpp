#include "numpy.h"

// include FFT library
#include <fftw3.h>

std::vector<std::vector<int>> Numpy::Combinations(int n, int k) {
    // create the bitmask that will get permuted and the resulting vector
    std::string bitmask(k, 1); bitmask.resize(n, 0); std::vector<std::vector<int>> combs;
 
    // generate the combinations
    do {std::vector<int> comb; comb.reserve(k);
        for (int j = 0; j < n; j++) {
            if (bitmask[j]) comb.push_back(j);
        } combs.push_back(comb);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    // return the result
    return combs;
}

Matrix<> Numpy::Repeat(const Matrix<>& A, int count, int axis) {
    // create the result matrix
    Matrix<> B(axis ? A.rows() : count * A.rows(), axis ? count * A.cols() : A.cols());

    // throw error if axis is out of bounds
    if (axis != 0 && axis != 1) throw std::runtime_error("UNKNOWN AXIS IN MATRIX REPEAT");

    // repeat rows
    if (axis == 0) {
        for (int i = 0; i < A.rows(); i++) {
            for (int j = 0; j < count; j++) {
                B.row(i * count + j) = A.row(i);
            }
        }
    }

    // repeat columns
    else if (axis == 1) {
        for (int i = 0; i < A.cols(); i++) {
            for (int j = 0; j < count; j++) {
                B.col(i * count + j) = A.col(i);
            }
        }
    }

    // return result
    return B;
}

Tensor<> Numpy::Kron(const Matrix<>& A, const Tensor<>& B) {
    // define the resulting tensor
    Tensor<> C(B.dimension(0), B.dimension(1), A.rows() * B.dimension(2), A.cols() * B.dimension(3));

    // perform the Kronecker product
    for (int i = 0; i < C.dimension(0); i++) {
        for (int j = 0; j < C.dimension(1); j++) {
            for (int k = 0; k < C.dimension(2); k++) {
                for (int l = 0; l < C.dimension(3); l++) {
                    C(i, j, k, l) = A(k / B.dimension(2), l / B.dimension(3)) * B(i, j, k % B.dimension(2), l % B.dimension(3));
                }
            }
        }
    }

    // return the result
    return C;
}

Vector<std::complex<double>> Numpy::FFT(Vector<std::complex<double>> in, int sign) {
    // define the output vector
    Vector<std::complex<double>> out(in.size());

    // plan the FFT
    fftw_plan plan = fftw_plan_dft_1d(in.size(), reinterpret_cast<fftw_complex*>(in.data()), reinterpret_cast<fftw_complex*>(out.data()), sign, FFTW_ESTIMATE);

    // execute the FFT and destroy the plan
    fftw_execute(plan); fftw_destroy_plan(plan);

    // return output
    return out / (sign == 1 ? in.size() : 1);
}

Vector<std::complex<double>> Numpy::FFT(Vector<std::complex<double>> in) {
    return FFT(in, FFTW_FORWARD);
}

Vector<std::complex<double>> Numpy::IFFT(Vector<std::complex<double>> in) {
    return FFT(in, FFTW_BACKWARD);
}

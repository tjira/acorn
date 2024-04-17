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

void Numpy::HFFT(std::complex<double>* in, double* out, const std::vector<int>& sizes) {
    fftw_plan plan = fftw_plan_dft_c2r(sizes.size(), sizes.data(), reinterpret_cast<fftw_complex*>(in), out, FFTW_ESTIMATE); fftw_execute(plan); fftw_destroy_plan(plan);
}

void Numpy::FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes, int sign) {
    fftw_plan plan = fftw_plan_dft(sizes.size(), sizes.data(), reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), sign, FFTW_ESTIMATE); fftw_execute(plan); fftw_destroy_plan(plan);
}

Matrix<double> Numpy::HFFT(Matrix<std::complex<double>> in) {
    // extract dimension
    int m = in.rows(), n = in.cols();

    // initialize the output vector and perform the FFT
    std::vector<double> out((in.size() - 1) * 2); Numpy::HFFT(in.data(), out.data(), std::vector<int>{m, n});

    // transform the output to the correct format and return
    return Eigen::Map<Vector<double>>(out.data(), in.size());
}

Matrix<std::complex<double>> Numpy::FFT(Matrix<std::complex<double>> in, int sign) {
    // extract dimension
    int m = in.rows(), n = in.cols();

    // initialize the output vector and perform the FFT
    std::vector<std::complex<double>> out(in.size()); Numpy::FFT(in.data(), out.data(), std::vector<int>{m, n}, sign);

    // transform the output to the correct format and return
    return Eigen::Map<Matrix<std::complex<double>>>(out.data(), in.rows(), in.cols()) / (sign == 1 ? out.size() : 1);
}

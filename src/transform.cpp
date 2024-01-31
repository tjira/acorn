#include "transform.h"

Matrix<> Transform::Single(const Matrix<>& A, const Matrix<>& C) {
    // create the transformed matrix
    Matrix<> Amo(A.rows(), A.cols());

    // perform the transformation
    #pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            for (int a = 0; a < Amo.rows(); a++) {
                for (int b = 0; b < Amo.cols(); b++) {
                    Amo(i, j) += A(a, b) * C(a, i) * C(b, j);
                }
            }
        }
    }

    // return the transformed matrix
    return Amo;
}

Matrix<> Transform::SingleSpin(const Matrix<>& A, const Matrix<>& C) {
    // expand the dimentsions of the MO basis matrix
    Matrix<> Ams = Numpy::Repeat(Numpy::Repeat(Transform::Single(A, C), 2, 0), 2, 1);

    // create repeating zeros and ones and tile it to a matrix
    Vector<double> ind(Ams.rows()); std::iota(ind.begin(), ind.end(), 0); ind = ind.unaryExpr([](auto s) {return double(int(s) % 2);});
    Matrix<bool> indm(ind.size(), ind.size()); for (size_t i = 0; i < ind.size(); i++) indm.row(i) = ind.transpose().array() == ind(i);

    // element wise multiply and return
    return Ams.array() * indm.cast<double>().array();
}

Tensor<> Transform::Coulomb(const Tensor<>& J, const Matrix<>& C) {
    // declare the Coulomb integral in molecular orbital basis and tensors of partial transform
    Tensor<4> J01(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J01.setZero();
    Tensor<4> J02(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J02.setZero();
    Tensor<4> J03(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J03.setZero();
    Tensor<4> Jmo(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); Jmo.setZero();

    // first 5th order partial transform
    #pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < J.dimension(0); i++) {
        for (int a = 0; a < J.dimension(0); a++) {
            for (int b = 0; b < J.dimension(1); b++) {
                for (int c = 0; c < J.dimension(2); c++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        J01(i, b, c, d) += J(a, b, c, d) * C(a, i);
                    }
                }
            }
        }
    }

    // second 5th order partial transform
    #pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < J.dimension(0); i++) {
        for (int j = 0; j < J.dimension(1); j++) {
            for (int b = 0; b < J.dimension(1); b++) {
                for (int c = 0; c < J.dimension(2); c++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        J02(i, j, c, d) += J01(i, b, c, d) * C(b, j);
                    }
                }
            }
        }
    }

    // third 5th order partial transform
    #pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < J.dimension(0); i++) {
        for (int j = 0; j < J.dimension(1); j++) {
            for (int k = 0; k < J.dimension(2); k++) {
                for (int c = 0; c < J.dimension(2); c++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        J03(i, j, k, d) += J02(i, j, c, d) * C(c, k);
                    }
                }
            }
        }
    }

    // fourth 5th order partial transform
    #pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < J.dimension(0); i++) {
        for (int j = 0; j < J.dimension(1); j++) {
            for (int k = 0; k < J.dimension(2); k++) {
                for (int l = 0; l < J.dimension(3); l++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        Jmo(i, j, k, l) += J03(i, j, k, d) * C(d, l);
                    }
                }
            }
        }
    }

    // return the integral in molecular orbital basis
    return Jmo;
}

Tensor<> Transform::CoulombSpin(const Tensor<>& J, const Matrix<>& C) {
    // calculate the integrals in the MO basis
    Matrix<> Jmo = Eigen::Map<Matrix<>>(Coulomb(J, C).data(), C.rows() * C.rows(), C.rows() * C.rows());

    // create the MS masis integrals container
    Matrix<> Jms(4 * C.cols() * C.cols(), 4 * C.rows() * C.rows()); Jms.setZero();

    // perform the transform
    #pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < 2 * C.rows(); i++) {
        for (int j = 0; j < 2 * C.rows(); j++) {
            for (int k = 0; k < 2 * C.rows(); k++) {
                for (int l = 0; l < 2 * C.rows(); l++) {
                    Jms(i * 2 * C.rows() + j, k * 2 * C.rows() + l) = Jmo(i / 2 * C.rows() + k / 2, j / 2 * C.rows() + l / 2) * (i % 2 == k % 2) * (j % 2 == l % 2);
                }
            }
        }
    }

    // return the resulting integrals in the MS basis
    return Eigen::TensorMap<Tensor<4>>(Jms.data(), 2 * C.rows(), 2 * C.rows(), 2 * C.rows(), 2 * C.rows()).shuffle(Eigen::array<int, 4>{3, 2, 1, 0});
}

#include "transform.h"

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

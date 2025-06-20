#include <unsupported/Eigen/MatrixFunctions>

#include <omp.h>
extern "C" {
    using namespace Eigen; typedef unsigned long ulong;

    void logm(double *B, double *A, const ulong dim) {
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AT(A, dim, dim);
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> BT(B, dim, dim);

        BT = AT.log();
    }
}

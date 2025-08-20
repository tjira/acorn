#include <unsupported/Eigen/MatrixFunctions>

extern "C" {
    using namespace Eigen; typedef unsigned long ulong;

    double det(double *A, const ulong dim) {
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AM(A, dim, dim);

        return AM.determinant();
    }

    void logm(double *B, double *A, const ulong dim) {
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AM(A, dim, dim);
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> BM(B, dim, dim);

        BM = AM.log();
    }
}

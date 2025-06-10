#define EIGEN_USE_THREADS

#include <unsupported/Eigen/MatrixFunctions>
#include    <unsupported/Eigen/CXX11/Tensor>

template<int RA, int RB, int CS>
void contract(double *C, const unsigned long *dC, double *A, const unsigned long *dA, double *B, const unsigned long *dB, const int *pairs) {
    Eigen::ThreadPool pool(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1);

    const int RC = RA + RB - 2 * CS;

    Eigen::DSizes<Eigen::Index,RA> dimsA; for(int i = 0; i < RA; i++) dimsA[i] = static_cast<Eigen::Index>(dA[i]);
    Eigen::DSizes<Eigen::Index,RB> dimsB; for(int i = 0; i < RB; i++) dimsB[i] = static_cast<Eigen::Index>(dB[i]);
    Eigen::DSizes<Eigen::Index,RC> dimsC; for(int i = 0; i < RC; i++) dimsC[i] = static_cast<Eigen::Index>(dC[i]);

    Eigen::TensorMap<Eigen::Tensor<double, RA, Eigen::RowMajor>> AT(A, dimsA);
    Eigen::TensorMap<Eigen::Tensor<double, RB, Eigen::RowMajor>> BT(B, dimsB);
    Eigen::TensorMap<Eigen::Tensor<double, RC, Eigen::RowMajor>> CT(C, dimsC);

    Eigen::array<Eigen::IndexPair<int>, CS> axes;

    for (int i = 0; i < CS; ++i) {
        axes[i] = Eigen::IndexPair<int>{pairs[2 * i], pairs[2 * i + 1]};
    }

    CT.device(Eigen::ThreadPoolDevice(&pool, pool.NumThreads())) = AT.contract(BT, axes);
}

#include <omp.h>
extern "C" {
    using namespace Eigen; typedef unsigned long ulong;

    void contract(double *C, const ulong *dC, ulong rC, double *A, const ulong *dA, ulong rA, double *B, const ulong *dB, ulong rB, const int *pairs, ulong npairs) {
        if (rA == 2 && rB == 2 && npairs == 1) return contract<2, 2, 1>(C, dC, A, dA, B, dB, pairs);
        if (rA == 2 && rB == 4 && npairs == 1) return contract<2, 4, 1>(C, dC, A, dA, B, dB, pairs);
        if (rA == 2 && rB == 4 && npairs == 2) return contract<2, 4, 2>(C, dC, A, dA, B, dB, pairs);

        throw std::invalid_argument("UNSUPPORTED CONTRACTION SIGNATURE");
    }

    void logm(double *B, double *A, const ulong n) {
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AT(A, n, n);
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> BT(B, n, n);

        BT = AT.log();
    }
}

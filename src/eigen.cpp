#define EIGEN_USE_THREADS

#include <unsupported/Eigen/CXX11/Tensor>

#include <omp.h>
extern "C" {
    using namespace Eigen;

    void contract_221(double *C, const unsigned long *dC, double *A, const unsigned long *dA, double *B, const unsigned long *dB, const int *pairs) {
        ThreadPool pool(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1);

        TensorMap<Tensor<double, 2, RowMajor>> AT(A, (long)dA[0], (long)dA[1]);
        TensorMap<Tensor<double, 2, RowMajor>> BT(B, (long)dB[0], (long)dB[1]);
        TensorMap<Tensor<double, 2, RowMajor>> CT(C, (long)dC[0], (long)dC[1]);

        array<IndexPair<int>, 1> axes = {
            IndexPair<int>{pairs[0], pairs[1]}
        };

        CT.device(ThreadPoolDevice(&pool, pool.NumThreads())) = AT.contract(BT, axes);
    }

    void contract_241(double *C, const unsigned long *dC, double *A, const unsigned long *dA, double *B, const unsigned long *dB, const int *pairs) {
        ThreadPool pool(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1);

        TensorMap<Tensor<double, 2, RowMajor>> AT(A, (long)dA[0], (long)dA[1]                          );
        TensorMap<Tensor<double, 4, RowMajor>> BT(B, (long)dB[0], (long)dB[1], (long)dB[2], (long)dB[3]);
        TensorMap<Tensor<double, 4, RowMajor>> CT(C, (long)dC[0], (long)dC[1], (long)dC[2], (long)dC[3]);

        array<IndexPair<int>, 1> axes = {
            IndexPair<int>{pairs[0], pairs[1]}
        };

        CT.device(ThreadPoolDevice(&pool, pool.NumThreads())) = AT.contract(BT, axes);
    }

    void contract_242(double *C, const unsigned long *dC, double *A, const unsigned long *dA, double *B, const unsigned long *dB, const int *pairs) {
        ThreadPool pool(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1);

        TensorMap<Tensor<double, 2, RowMajor>> AT(A, (long)dA[0], (long)dA[1]                          );
        TensorMap<Tensor<double, 4, RowMajor>> BT(B, (long)dB[0], (long)dB[1], (long)dB[2], (long)dB[3]);
        TensorMap<Tensor<double, 2, RowMajor>> CT(C, (long)dC[0], (long)dC[1]                          );

        array<IndexPair<int>, 2> axes = {
            IndexPair<int>{pairs[0], pairs[1]},
            IndexPair<int>{pairs[2], pairs[3]}
        };

        CT.device(ThreadPoolDevice(&pool, pool.NumThreads())) = AT.contract(BT, axes);
    }
}

#define EIGEN_USE_THREADS

#include <boost/preprocessor/seq/for_each.hpp>
#include   <boost/preprocessor/tuple/elem.hpp>

#include <unsupported/Eigen/MatrixFunctions>
#include    <unsupported/Eigen/CXX11/Tensor>

#define CONTRACT_CASES ((2,2,1))((2,4,1))((2,4,2))

#define EMIT_CONTRACT_CASE(r, data, elem) \
    if (rA == BOOST_PP_TUPLE_ELEM(3, 0, elem) && rB == BOOST_PP_TUPLE_ELEM(3, 1, elem) && npairs == BOOST_PP_TUPLE_ELEM(3, 2, elem)) { \
        return contract<BOOST_PP_TUPLE_ELEM(3, 0, elem), BOOST_PP_TUPLE_ELEM(3, 1, elem), BOOST_PP_TUPLE_ELEM(3, 2, elem)>(C, dC, A, dA, B, dB, pairs); \
    }

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
        BOOST_PP_SEQ_FOR_EACH(EMIT_CONTRACT_CASE, _, CONTRACT_CASES) throw std::invalid_argument("UNSUPPORTED CONTRACTION SIGNATURE");
    }

    void logm(double *B, double *A, const ulong n) {
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AT(A, n, n);
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> BT(B, n, n);

        BT = AT.log();
    }
}

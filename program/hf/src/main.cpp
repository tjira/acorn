#include "argparse.hpp"
#include "tensor.h"

inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> toMatrix(Eigen::Tensor<double, 2, Eigen::ColMajor> A) {return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(A.data(), A.dimension(0), A.dimension(1));}
inline Eigen::Tensor<double, 2, Eigen::ColMajor> toTensor(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> A) {return Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(A.data(), A.rows(), A.cols());}

#include <unsupported/Eigen/MatrixFunctions>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none);

    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    Matrix<> V = Matrix<>::Load("V.mat"), T = Matrix<>::Load("T.mat"), S = Matrix<>::Load("S.mat"); Tensor<> J = Tensor<>::Load("J.mat");

    double E = 0; int nocc = 5;

    Tensor<> Janti = J - J.transpose({0, 3, 2, 1}) * 0.5;

    // Matrix<> H = T + V, F = H, C(S.rows(), S.cols()), D(S.rows(), S.cols());

    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> H = T.mat + V.mat, F = H, C(S.rows(), S.cols()), D(S.rows(), S.cols()), e;
    Eigen::Tensor<double, 4, Eigen::ColMajor> ERI = Janti.ten;

    for (int i = 0; i < 100; i++) {
        // calculate the Fock matrix
        F = H + toMatrix(ERI.contract(toTensor(D), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> solver(F, S.mat);

        // exteract the eigenvalues and eigenvectors
        C = solver.eigenvectors(), e = solver.eigenvalues();

        // save previous density matrix and energy
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Dp = D; double Ep = E;

        // calculate the new density matrix
        D = 2 * C.leftCols(nocc) * C.leftCols(nocc).transpose();

        // calculate the new energy
        E = 0.5 * D.cwiseProduct(H + F).sum();

        // calculate the E and D errors
        double Eerr = std::abs(E - Ep), Derr = (D - Dp).norm();

        // finish if covergence reached
        if (Eerr < 1e-8 && Derr < 1e-8) {break;}

        std::cout << E << std::endl;
    }

    std::cout << "Energy: " << E << std::endl;

}

#include "argparse.hpp"

#include "system.h"
#include "tensor.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Number of SCF iterations.").default_value(100).scan<'i', int>();
    program.add_argument("-t", "--threshold").help("-- Convergence threshold.").default_value(1e-8).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // load the system and integrals from disk
    System system("molecule.xyz"); Matrix<> V = Matrix<>::Load("V.mat"), T = Matrix<>::Load("T.mat"), S = Matrix<>::Load("S.mat"); Tensor<> J = Tensor<>::Load("J.mat");

    // initialize all the matrices used throughout the SCF procedure and the energy
    Matrix<> H = T + V, F = H, C(S.rows(), S.cols()), D(S.rows(), S.cols()), eps(S.rows(), 1);

    // initialize the coulomb contraction axes, the energy placeholder and the number of occupied orbitals
    Eigen::IndexPair<int> first(2, 0), second(3, 1); double E = 0; int nocc = system.nocc();

    // start the SCF procedure
    for (int i = 0; i < program.get<int>("-i"); i++) {
        // calculate the iteration Fock matrix
        F = H + (J - 0.5 * J.t({0, 3, 2, 1})).contract(D.tensor(), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}).matrix();

        // solve the eigenproblem and save the previous D and E values
        std::tie(eps, C) = F.eigh(S); Matrix Dp = D; double Ep = E;

        // calculate the new density matrix and energy
        D = 2 * C.leftcols(nocc).dot(C.leftcols(nocc).t()); E = 0.5 * (D * (H + F)).sum();

        // print the iteration info
        std::printf("%4d %20.14f %.2e %.2e %s\n", i + 1, E, std::abs(E - Ep), (D - Dp).norm(), "");

        // finish if covergence reached
        if (double thresh = program.get<double>("-t"); std::abs(E - Ep) < thresh && (D - Dp).norm() < thresh) break;
        else if (i == program.get<int>("-i") - 1) {
            throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE SCF REACHED.");
        }
    }

    // save the final matrices
    C.save("C.mat"), D.save("D.mat"), eps.save("eps.mat");

    // print the final energy
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n", E + system.nuclearRepulsion());
}

#include "system.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Number of SCF iterations.").default_value(100).scan<'i', int>();
    program.add_argument("-t", "--threshold").help("-- Convergence threshold.").default_value(1e-8).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the system from disk
    System system("molecule.xyz");

    // load the integrals in AO basis from disk
    tp = Timer::Now(); std::cout << "V_AO MATRIX READING: " << std::flush; EigenMatrix<> V = Eigen::LoadMatrix("V_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "T_AO MATRIX READING: " << std::flush; EigenMatrix<> T = Eigen::LoadMatrix("T_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "S_AO MATRIX READING: " << std::flush; EigenMatrix<> S = Eigen::LoadMatrix("S_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "J_AO TENSOR READING: " << std::flush; EigenTensor<> J = Eigen::LoadTensor("J_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // initialize all the matrices used throughout the SCF procedure and the energy
    EigenMatrix<> H = T + V, F = H, C(S.rows(), S.cols()), D(S.rows(), S.cols()), eps(S.rows(), 1);

    // initialize the contraction axes and the energy placeholder
    Eigen::IndexPair<int> first(2, 0), second(3, 1); double E = 0;

    // print the header
    std::printf("\n%6s %20s %8s %8s %12s\n", "ITER", "ENERGY", "|dE|", "|dD|", "TIME");

    // start the SCF procedure
    for (int i = 0; i < program.get<int>("-i"); i++) {

        // reset the timer
        tp = Timer::Now();

        // calculate the electron-electron repulsion
        EigenTensor<2> VEE = (J - 0.5 * J.shuffle(Eigen::array<int, 4>{0, 3, 2, 1})).contract(TENSORMAP(D), Eigen::array<Eigen::IndexPair<int>, 2>{first, second});

        // calculate the Fock matrix
        F = H + MATRIXMAP(VEE);

        // solve the eigenproblem and save the previous D and E values
        std::tie(eps, C) = Eigen::Gsaes(F, S); EigenMatrix<> Dp = D; double Ep = E;

        // calculate the new density
        D = 2 * C.leftCols(system.nocc()) * C.leftCols(system.nocc()).transpose();

        // calculate the new energy
        E = 0.5 * D.cwiseProduct(H + F).sum();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %.2e %s\n", i + 1, E, std::abs(E - Ep), (D - Dp).norm(), Timer::Format(Timer::Elapsed(tp)).c_str());

        // finish if covergence reached
        if (double thresh = program.get<double>("-t"); std::abs(E - Ep) < thresh && (D - Dp).norm() < thresh) {std::cout << std::endl; break;}
        else if (i == program.get<int>("-i") - 1) {
            throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE SCF REACHED.");
        }
    }

    // save the final matrices
    tp = Timer::Now(); std::cout << "E_MO VECTOR WRITING: " << std::flush; Eigen::Write("E_MO.mat", eps); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "C_MO MATRIX WRITING: " << std::flush; Eigen::Write("C_MO.mat", C); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "D_MO MATRIX WRITING: " << std::flush; Eigen::Write("D_MO.mat", D); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // print the final energy and total time
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", E + system.nuclearRepulsion(), Timer::Format(Timer::Elapsed(start)).c_str());
}

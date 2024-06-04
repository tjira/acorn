#include "system.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-d", "--diis").help("-- Size of the DIIS subspace.").default_value(0).scan<'i', int>();
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-i", "--iterations").help("-- Number of SCF iterations.").default_value(100).scan<'i', int>();
    program.add_argument("-t", "--threshold").help("-- Convergence threshold.").default_value(1e-8).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters
    int diis = program.get<int>("-d"), iters = program.get<int>("-i"); double thresh = program.get<double>("-t");

    // load the integrals in AO basis and system from disk
    MEASURE("SYSTEM AND INTEGRALS IN AO BASIS READING: ",
        Matrix    V = Eigen::LoadMatrix("V_AO.mat");
        Matrix    T = Eigen::LoadMatrix("T_AO.mat");
        Matrix    S = Eigen::LoadMatrix("S_AO.mat");
        Tensor<4> J = Eigen::LoadTensor("J_AO.mat");
        System system(program.get("-f"));
    )

    // initialize all the matrices used throughout the SCF procedure and the energy
    Matrix H = T + V, F = H, Cmo(S.rows(), S.cols()), Dmo(S.rows(), S.cols()), Emo(S.rows(), 1);

    // initialize the contraction axes, the energy placeholder and the Fock and error vector container
    Eigen::IndexPair<int> first(2, 0), second(3, 1); double Eel = 0; std::vector<Matrix> es, fs;

    // print the header
    std::printf("\n%6s %20s %8s %8s %12s\n", "ITER", "ENERGY", "|dE|", "|dD|", "TIME");

    // start the SCF procedure
    for (int i = 0; i < iters; i++) {

        // reset the timer
        tp = Timer::Now();

        // calculate the electron-electron repulsion
        Tensor<2> VEE = (J - 0.5 * J.shuffle(Eigen::array<int, 4>{0, 3, 2, 1})).contract(TENSORMAP(Dmo), Eigen::array<Eigen::IndexPair<int>, 2>{first, second});

        // calculate the Fock matrix and define previous values
        F = H + MATRIXMAP(VEE); Matrix Dmop = Dmo; double Eelp = Eel;

        // calculate the error vector and append it to the container along with the Fock matrix
        Matrix e = S * Dmo * F - F * Dmo * S; if (i) es.push_back(e), fs.push_back(F);

        // truncate the error and Fock vector containers
        if (i > diis) {es.erase(es.begin()), fs.erase(fs.begin());}

        // perform DIIS extrapolation
        if (diis && i >= diis) {

            // define the DIIS subspace matrices
            Matrix B = Matrix::Ones(diis + 1, diis + 1), b = Vector::Zero(diis + 1); B(diis, diis) = 0, b(diis) = 1;

            // fill the DIIS matrix
            for (int j = 0; j < diis; j++) for (int k = j; k < diis; k++) B(j, k) = es.at(j).cwiseProduct(es.at(k)).sum(), B(k, j) = B(j, k);

            // solve the DIIS equations
            Vector c = B.colPivHouseholderQr().solve(b);

            // extrapolate the Fock matrix
            F = c(0) * fs.at(0); for (int j = 1; j < diis; j++) F += c(j) * fs.at(j);
        }

        // solve the Roothaan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, S); Emo = solver.eigenvalues(), Cmo = solver.eigenvectors();

        // calculate the new density
        Dmo = 2 * Cmo.leftCols(system.nocc()) * Cmo.leftCols(system.nocc()).transpose();

        // calculate the new energy
        Eel = 0.5 * Dmo.cwiseProduct(H + F).sum();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %.2e %s %s\n", i + 1, Eel, std::abs(Eel - Eelp), (Dmo - Dmop).norm(), Timer::Format(Timer::Elapsed(tp)).c_str(), diis && i >= diis ? "DIIS" : "");

        // finish if covergence reached
        if (std::abs(Eel - Eelp) < thresh && (Dmo - Dmop).norm() < thresh) {std::cout << std::endl; break;}
        else if (i == iters - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE SCF REACHED.");
    }

    // save the final matrices
    MEASURE("HARTREE-FOCK RESULTS WRITING: ",
        Eigen::Write("C_MO.mat", Cmo);
        Eigen::Write("E_MO.mat", Emo);
        Eigen::Write("D_MO.mat", Dmo);
    )

    // print the final energy and total time
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", Eel + system.nuclearRepulsion(), Timer::Format(Timer::Elapsed(start)).c_str());
}

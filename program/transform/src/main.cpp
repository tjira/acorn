#include "transform.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Integral Transform Engine", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // load the integrals in AO basis from disk
    EigenMatrix<> V = Eigen::LoadMatrix("V_AO.mat");
    EigenMatrix<> T = Eigen::LoadMatrix("T_AO.mat");
    EigenMatrix<> S = Eigen::LoadMatrix("S_AO.mat");
    EigenMatrix<> C = Eigen::LoadMatrix("C_MO.mat");
    EigenTensor<> J = Eigen::LoadTensor("J_AO.mat");

    // load the molecular orbital energies from disk
    EigenMatrix<> eps = Eigen::LoadMatrix("E_MO.mat");

    // one-electron integrals in molecular spatial orbital basis
    Eigen::Write("V_MO.mat", Transform::SingleSpatial(V, C));
    Eigen::Write("T_MO.mat", Transform::SingleSpatial(T, C));
    Eigen::Write("S_MO.mat", Transform::SingleSpatial(S, C));

    // one-electron integrals in molecular spinorbital basis
    Eigen::Write("V_MS.mat", Transform::SingleSpin(V, C));
    Eigen::Write("T_MS.mat", Transform::SingleSpin(T, C));
    Eigen::Write("S_MS.mat", Transform::SingleSpin(S, C));

    // two-electron integrals in molecular spatial orbital basis
    Eigen::Write("J_MO.mat", Transform::CoulombSpatial(J, C));

    // two-electron integrals in molecular spinorbital basis
    Eigen::Write("J_MS.mat", Transform::CoulombSpin(J, C));
}

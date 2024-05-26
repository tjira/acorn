#include "transform.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Integral Transform Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the coefficient matrix
    MEASURE("COEFFICIENT MATRIX IN MO BASIS READING: ", EigenMatrix<> C = Eigen::LoadMatrix("C_MO.mat"))

    // print the new line
    std::cout << std::endl;

    // load the integrals in AO basis from disk
    MEASURE("NUCLEAR INTEGRALS IN AO BASIS READING: ", EigenMatrix<> V = Eigen::LoadMatrix("V_AO.mat"))
    MEASURE("KINETIC INTEGRALS IN AO BASIS READING: ", EigenMatrix<> T = Eigen::LoadMatrix("T_AO.mat"))
    MEASURE("OVERLAP INTEGRALS IN AO BASIS READING: ", EigenMatrix<> S = Eigen::LoadMatrix("S_AO.mat"))
    MEASURE("COULOMB INTEGRALS IN AO BASIS READING: ", EigenTensor<> J = Eigen::LoadTensor("J_AO.mat"))

    // print a new line
    std::cout << std::endl;

    // one-electron integrals in molecular spatial orbital basis calculation
    MEASURE("NUCLEAR INTEGRALS MO TRANSFORM: ", EigenMatrix<> Vmo = Transform::SingleSpatial(V, C))
    MEASURE("KINETIC INTEGRALS MO TRANSFORM: ", EigenMatrix<> Tmo = Transform::SingleSpatial(T, C))
    MEASURE("OVERLAP INTEGRALS MO TRANSFORM: ", EigenMatrix<> Smo = Transform::SingleSpatial(S, C))

    // one-electron integrals in molecular spinorbital basis calculation
    MEASURE("NUCLEAR INTEGRALS MS TRANSFORM: ", EigenMatrix<> Vms = Transform::SingleSpin(V, C))
    MEASURE("KINETIC INTEGRALS MS TRANSFORM: ", EigenMatrix<> Tms = Transform::SingleSpin(T, C))
    MEASURE("OVERLAP INTEGRALS MS TRANSFORM: ", EigenMatrix<> Sms = Transform::SingleSpin(S, C))

    // two-electron integrals in molecular spatial orbital basis calculation
    MEASURE("COULOMB INTEGRALS MO TRANSFORM: ", EigenTensor<> Jmo = Transform::CoulombSpatial(J, C))

    // two-electron integrals in molecular spinorbital basis calculation
    MEASURE("COULOMB INTEGRALS MS TRANSFORM: ", EigenTensor<> Jms = Transform::CoulombSpin(J, C))

    // print new line
    std::cout << std::endl;

    // one-electron integrals in molecular spatial orbital basis writing
    MEASURE("NUCLEAR INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("V_MO.mat", Vmo))
    MEASURE("KINETIC INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("T_MO.mat", Tmo))
    MEASURE("OVERLAP INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("S_MO.mat", Smo))

    // one-electron integrals in molecular spinorbital basis writing
    MEASURE("NUCLEAR INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("V_MS.mat", Vms))
    MEASURE("KINETIC INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("T_MS.mat", Tms))
    MEASURE("OVERLAP INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("S_MS.mat", Sms))

    // two-electron integrals in molecular spatial orbital basis writing
    MEASURE("COULOMB INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("J_MO.mat", Jmo))

    // two-electron integrals in molecular spinorbital basis writing
    MEASURE("COULOMB INTEGRALS IN MS BASIS WRITING: ", Eigen::Write("J_MS.mat", Jms))

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

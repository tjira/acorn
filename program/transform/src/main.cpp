#include "transform.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Integral Transform Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--orbital").help("-- Transforms to the basis of spatial molecular orbitals.").default_value(false).implicit_value(true);
    program.add_argument("-s", "--spinorbital").help("-- Transforms to the basis of molecular spinorbitals.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // define all the matrices and tensors used throughout the program
    EigenMatrix<> C, V, T, S, Vmo, Tmo, Smo, Vms, Tms, Sms; EigenTensor<> J, Jmo, Jms;

    // load the coefficient matrix
    if (program.get<bool>("-o") || program.get<bool>("-s")) {MEASURE("COEFFICIENT MATRIX IN MO BASIS READING: ", C = Eigen::LoadMatrix("C_MO.mat"))}

    // print the new line
    if (program.get<bool>("-o") || program.get<bool>("-s")) std::cout << std::endl;

    // load the integrals in AO basis from disk
    if (program.get<bool>("-o") || program.get<bool>("-s")) {MEASURE("NUCLEAR INTEGRALS IN AO BASIS READING: ", V = Eigen::LoadMatrix("V_AO.mat"))}
    if (program.get<bool>("-o") || program.get<bool>("-s")) {MEASURE("KINETIC INTEGRALS IN AO BASIS READING: ", T = Eigen::LoadMatrix("T_AO.mat"))}
    if (program.get<bool>("-o") || program.get<bool>("-s")) {MEASURE("OVERLAP INTEGRALS IN AO BASIS READING: ", S = Eigen::LoadMatrix("S_AO.mat"))}
    if (program.get<bool>("-o") || program.get<bool>("-s")) {MEASURE("COULOMB INTEGRALS IN AO BASIS READING: ", J = Eigen::LoadTensor("J_AO.mat"))}

    // print a new line
    if (program.get<bool>("-o") || program.get<bool>("-s")) std::cout << std::endl;

    // one-electron integrals in molecular spatial orbital basis calculation
    if (program.get<bool>("-o")) {MEASURE("NUCLEAR INTEGRALS MO TRANSFORM: ", Vmo = Transform::SingleSpatial(V, C))}
    if (program.get<bool>("-o")) {MEASURE("KINETIC INTEGRALS MO TRANSFORM: ", Tmo = Transform::SingleSpatial(T, C))}
    if (program.get<bool>("-o")) {MEASURE("OVERLAP INTEGRALS MO TRANSFORM: ", Smo = Transform::SingleSpatial(S, C))}

    // one-electron integrals in molecular spinorbital basis calculation
    if (program.get<bool>("-s")) {MEASURE("NUCLEAR INTEGRALS MS TRANSFORM: ", Vms = Transform::SingleSpin(V, C))}
    if (program.get<bool>("-s")) {MEASURE("KINETIC INTEGRALS MS TRANSFORM: ", Tms = Transform::SingleSpin(T, C))}
    if (program.get<bool>("-s")) {MEASURE("OVERLAP INTEGRALS MS TRANSFORM: ", Sms = Transform::SingleSpin(S, C))}

    // two-electron integrals in molecular spatial orbital basis calculation
    if (program.get<bool>("-o")) {MEASURE("COULOMB INTEGRALS MO TRANSFORM: ", Jmo = Transform::CoulombSpatial(J, C))}

    // two-electron integrals in molecular spinorbital basis calculation
    if (program.get<bool>("-s")) {MEASURE("COULOMB INTEGRALS MS TRANSFORM: ", Jms = Transform::CoulombSpin(J, C))}

    // print new line
    if (program.get<bool>("-o") || program.get<bool>("-s")) std::cout << std::endl;

    // one-electron integrals in molecular spatial orbital basis writing
    if (program.get<bool>("-o")) {MEASURE("NUCLEAR INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("V_MO.mat", Vmo))}
    if (program.get<bool>("-o")) {MEASURE("KINETIC INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("T_MO.mat", Tmo))}
    if (program.get<bool>("-o")) {MEASURE("OVERLAP INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("S_MO.mat", Smo))}

    // one-electron integrals in molecular spinorbital basis writing
    if (program.get<bool>("-s")) {MEASURE("NUCLEAR INTEGRALS IN MS BASIS WRITING: ", Eigen::Write("V_MS.mat", Vms))}
    if (program.get<bool>("-s")) {MEASURE("KINETIC INTEGRALS IN MS BASIS WRITING: ", Eigen::Write("T_MS.mat", Tms))}
    if (program.get<bool>("-s")) {MEASURE("OVERLAP INTEGRALS IN MS BASIS WRITING: ", Eigen::Write("S_MS.mat", Sms))}

    // two-electron integrals in molecular spatial orbital basis writing
    if (program.get<bool>("-o")) {MEASURE("COULOMB INTEGRALS IN MO BASIS WRITING: ", Eigen::Write("J_MO.mat", Jmo))}

    // two-electron integrals in molecular spinorbital basis writing
    if (program.get<bool>("-s")) {MEASURE("COULOMB INTEGRALS IN MS BASIS WRITING: ", Eigen::Write("J_MS.mat", Jms))}

    // print new line
    if (program.get<bool>("-s")) std::cout << std::endl;

    // orbital energies in molecular spinorbital basis calculation
    if (program.get<bool>("-s")) {MEASURE("ORBITAL ENERGIES READ, MS TRANSFORM AND WRITE: ", Eigen::Write("E_MS.mat", Eigen::Repeat(Eigen::LoadMatrix("E_MO.mat"), 2, 0)))}

    // print the total time
    std::cout << (program.get<bool>("-o") || program.get<bool>("-s") ? "\n" : "") << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

#include "integral.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-b", "--basis").help("-- Basis set used for integrals.").default_value("STO-3G");
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the system from disk and initialize the basis set
    std::ifstream fstream(program.get("-f")); if (!fstream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");
    std::vector<libint2::Atom> atoms = libint2::read_dotxyz(fstream); libint2::BasisSet shells(program.get("-b"), atoms);

    // print the timing of the system initialization
    std::cout << "SYSTEM INITIALIZATION: " << Timer::Format(Timer::Elapsed(tp)) << std::endl << std::endl;

    // calculate the integrals
    libint2::initialize();
    MEASURE("NUCLEAR INTEGRALS IN AO BASIS CALCULATION: ", EigenMatrix<> V = Integral::Nuclear(atoms, shells))
    MEASURE("KINETIC INTEGRALS IN AO BASIS CALCULATION: ", EigenMatrix<> T = Integral::Kinetic(shells))
    MEASURE("OVERLAP INTEGRALS IN AO BASIS CALCULATION: ", EigenMatrix<> S = Integral::Overlap(shells))
    MEASURE("COULOMB INTEGRALS IN AO BASIS CALCULATION: ", EigenTensor<> J = Integral::Coulomb(shells))
    libint2::finalize();

    // new line
    std::cout << std::endl;
    
    // save the integrals to disk
    MEASURE("NUCLEAR INTEGRALS IN AO BASIS WRITING: ", Eigen::Write("V_AO.mat", V))
    MEASURE("KINETIC INTEGRALS IN AO BASIS WRITING: ", Eigen::Write("T_AO.mat", T))
    MEASURE("OVERLAP INTEGRALS IN AO BASIS WRITING: ", Eigen::Write("S_AO.mat", S))
    MEASURE("COULOMB INTEGRALS IN AO BASIS WRITING: ", Eigen::Write("J_AO.mat", J))

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

#include "integral.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-b", "--basis").help("-- Basis set used for integrals.").default_value("STO-3G");
    program.add_argument("-f", "--file").help("-- System file in .xyz format.").default_value("molecule.xyz");
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
    tp = Timer::Now(); std::cout << "V_AO MATRIX CALCULATION: " << std::flush; EigenMatrix<> V = Integral::Nuclear(atoms, shells); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "T_AO MATRIX CALCULATION: " << std::flush; EigenMatrix<> T = Integral::Kinetic(shells); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "S_AO MATRIX CALCULATION: " << std::flush; EigenMatrix<> S = Integral::Overlap(shells); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "J_AO TENSOR CALCULATION: " << std::flush; EigenTensor<> J = Integral::Coulomb(shells); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    libint2::finalize();

    // new line
    std::cout << std::endl;
    
    // save the integrals to disk
    tp = Timer::Now(); std::cout << "V_AO MATRIX WRITING: " << std::flush; Eigen::Write("V_AO.mat", V); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "T_AO MATRIX WRITING: " << std::flush; Eigen::Write("T_AO.mat", T); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "S_AO MATRIX WRITING: " << std::flush; Eigen::Write("S_AO.mat", S); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "J_AO TENSOR WRITING: " << std::flush; Eigen::Write("J_AO.mat", J); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

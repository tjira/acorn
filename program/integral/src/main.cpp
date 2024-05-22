#include "argparse.hpp"
#include "integral.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Engine", "1.0", argparse::default_arguments::none);

    program.add_argument("-b", "--basis").help("-- Basis set used for integrals.").default_value("STO-3G");
    program.add_argument("-f", "--file").help("-- System file in .xyz format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    std::ifstream fstream(program.get("-f")); if (!fstream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");
    std::vector<libint2::Atom> atoms = libint2::read_dotxyz(fstream); libint2::BasisSet shells(program.get("-b"), atoms);

    libint2::initialize();
    Matrix<> V = Integral::Nuclear(atoms, shells); V.save("V.mat");
    Matrix<> T = Integral::Kinetic(shells); T.save("T.mat");
    Matrix<> S = Integral::Overlap(shells); S.save("S.mat");
    Tensor<> J = Integral::Coulomb(shells); J.save("J.mat");
    libint2::finalize();
}

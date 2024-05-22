#include "argparse.hpp"

#include "transform.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Integral Transform Engine", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // load the integrals in atomic spatial orbital basis from disk
    Matrix<> V = Matrix<>::Load("V.mat"), T = Matrix<>::Load("T.mat"), S = Matrix<>::Load("S.mat"), C = Matrix<>::Load("C.mat"); Tensor<> J = Tensor<>::Load("J.mat");

    Transform::Coulomb(J, C).save("Jmo.mat");
}

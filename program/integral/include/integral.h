#pragma once

#include "tensor.h"
#include <argparse.hpp>
#include <libint2.hpp>
#include <filesystem>

namespace Acorn {
    namespace Integral {
        // options struct
        struct Options {
            std::string basis, executable, file;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);

        // general integrals of the form (i|O|j) and (ij|O|kl) in chemists' notation
        torch::Tensor Single(libint2::Engine& engine, const libint2::BasisSet& shells);
        torch::Tensor Double(libint2::Engine& engine, const libint2::BasisSet& shells);

        // nuclear integral with the operator dependent of system coordinates
        torch::Tensor Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells);

        // kinetic, overlap, and coulomb integrals
        torch::Tensor Kinetic(const libint2::BasisSet& shells);
        torch::Tensor Overlap(const libint2::BasisSet& shells);
        torch::Tensor Coulomb(const libint2::BasisSet& shells);
    }
}

inline std::tuple<Acorn::Integral::Options, std::vector<timepoint>> Acorn::Integral::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Engine", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<timepoint> timers(2, std::chrono::high_resolution_clock().now()); auto tp = timers.at(0);

    // add the command line arguments
    program.add_argument("-b", "--basis").help("-- Basis set used for integrals.").default_value("STO-3G");
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // initialize the options
    Options opt = {program.get<std::string>("-b"), argv[0], program.get<std::string>("-f")};

    // return the options and timers
    return {opt, timers};
}

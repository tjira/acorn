#pragma once

#include "system.h"
#include <argparse.hpp>

namespace Acorn {
    namespace HF {
        // option struct
        struct Options {
            int diis, iters; double thresh; System system;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);
    }
}

inline std::tuple<Acorn::HF::Options, std::vector<timepoint>> Acorn::HF::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-d", "--diis").help("-- Size of the DIIS subspace.").default_value(0).scan<'i', int>();
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-i", "--iterations").help("-- Number of SCF iterations.").default_value(100).scan<'i', int>();
    program.add_argument("-t", "--threshold").help("-- Convergence threshold.").default_value(1e-8).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    Acorn::HF::Options opt = {program.get<int>("-d"), program.get<int>("-i"), program.get<double>("-t"), System(program.get("-f"))};

    // return the options and timers
    return {opt, timers};
}

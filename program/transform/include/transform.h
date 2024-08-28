#pragma once

#include "linalg.h"
#include "tensor.h"
#include <argparse.hpp>

namespace Acorn {
    namespace Transform {
        // options struct
        struct Options {
            bool orbital, spinorbital;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);

        // two electron integrals
        Tensor<4> CoulombSpin(const Tensor<4>& Jao, const Matrix& Cmo);
        Tensor<4> CoulombSpatial(const Tensor<4>& Jao, Matrix& Cmo);

        // one electron integrals
        Matrix SingleSpin(const Matrix& Aao, const Matrix& Cmo);
        Matrix SingleSpatial(const Matrix& Aao, Matrix& Cmo);
    }
}

inline std::tuple<Acorn::Transform::Options, std::vector<timepoint>> Acorn::Transform::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Transform Engine", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<timepoint> timers(2, std::chrono::high_resolution_clock().now()); auto tp = timers.at(0);

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--orbital").help("-- Transforms to the basis of spatial molecular orbitals.").default_value(false).implicit_value(true);
    program.add_argument("-s", "--spinorbital").help("-- Transforms to the basis of molecular spinorbitals.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // initialize the options
    Options opt = {program.get<bool>("-o"), program.get<bool>("-s")};

    // return the options and timers
    return {opt, timers};
}

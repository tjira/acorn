#pragma once

#include "tensor.h"
#include <argparse.hpp>

namespace Acorn {
    namespace MBPT {
        // options struct
        struct Options {
            int order;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);

        // custom MBPT functions
        double evaluate(const std::string& contrstr, const torch::Tensor& Jmsa, const torch::Tensor& Ems, int nos, int o);

        // number of threads
        inline int nthread = 1;
    }
}

inline std::tuple<Acorn::MBPT::Options, std::vector<timepoint>> Acorn::MBPT::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Moller-Plesset Perturbation Theory Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(3, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--order").help("-- Order of the theory.").default_value(2).scan<'i', int>();
    program.add_argument("-n", "--nthreads").help("-- Number of threads to use.").default_value(Acorn::MBPT::nthread).scan<'i', int>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set the number of threads
    Acorn::MBPT::nthread = program.get<int>("-n");

    // create the options
    Acorn::MBPT::Options opt = {program.get<int>("-o")};

    // return the options and timers
    return {opt, timers};
}

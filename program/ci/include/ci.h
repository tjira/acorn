#pragma once

#include "tensor.h"
#include <argparse.hpp>

namespace Acorn {
    namespace CI {
        // options struct
        struct Options {};

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);

        // number of threads
        inline int nthread = 1;
    }
}

inline std::tuple<Acorn::CI::Options, std::vector<timepoint>> Acorn::CI::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Configuration Interaction Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    //define the options
    Acorn::CI::Options opt = {};

    // return the options and timers
    return {opt, timers};
}

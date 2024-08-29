#pragma once

#include "expression.h"
#include "lz.h"
#include <argparse.hpp>
#include <random>

namespace Acorn {
    namespace CDYN {
        // options struct
        struct Options {
            std::vector<std::string> potential; double mass, momentum, position, step; int excstate, iters, log, nstate, seed, trajs; bool adiabatic, savetraj;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);

        // number of threads
        inline int nthread = 1;
    }
}

inline std::tuple<Acorn::CDYN::Options, std::vector<timepoint>> Acorn::CDYN::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Classical Dynamics Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-c", "--coordinate").help("-- Initial position of the system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-e", "--excstate").help("-- Initial diabatic state of the system.").default_value(0).scan<'i', int>();
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-l", "--log").help("-- Log interval for the simulation.").default_value(10).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-n", "--nthreads").help("-- Number of threads to use.").default_value(Acorn::CDYN::nthread).scan<'i', int>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-r", "--random").help("-- Random seed for the simulation.").default_value(1).scan<'i', int>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("-t", "--trajectories").help("-- Number of trajectories to run.").default_value(1).scan<'i', int>();
    program.add_argument("-u", "--potential").help("-- The potential curves in diabatic representation.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<std::string>{"0.5*x^2"});
    program.add_argument("--adiabatic").help("-- Perform the dynamics in adiabatic basis.").default_value(false).implicit_value(true);
    program.add_argument("--savetraj").help("-- Save the trajectories.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set the number of threads end extract the number of states
    Acorn::CDYN::nthread = program.get<int>("-n"); int nstate = std::sqrt(program.get<std::vector<std::string>>("-u").size());

    // extract the command line parameters
    Acorn::CDYN::Options opt = {
        program.get<std::vector<std::string>>("-u"), program.get<double>("-m"         ), program.get<double>("-p"        ), program.get<double>("-c"), program.get<double>("-s"),
        program.get<int                     >("-e"), program.get<int   >("-i"         ), program.get<int   >("-l"        ), nstate,                    program.get<int   >("-r"),
        program.get<int                     >("-t"), program.get<bool  >("--adiabatic"), program.get<bool  >("--savetraj")
    };

    // return the options and timers
    return {opt, timers};
}

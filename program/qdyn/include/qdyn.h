#pragma once

#include "fourier.h"
#include "tensor.h"
#include "wavefunction.h"
#include <argparse.hpp>

namespace Acorn {
    namespace QDYN {
        // options struct
        struct Options {
            double factor, mass, momentum, step; int dim, iters, optstates; bool adiabatic, align, imaginary, savewfn;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);
    }
}

inline std::tuple<Acorn::QDYN::Options, std::vector<timepoint>> Acorn::QDYN::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Quantum Dynamics Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-d", "--dimension").help("-- Dimension of the provided potential and wavefunction.").default_value(1).scan<'i', int>();
    program.add_argument("-f", "--factor").help("-- Factor for scaling the saved wavefunction.").default_value(1.0).scan<'g', double>();
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-o", "--optimize").help("-- Perform the dynamics in imaginary time for several orthogonal states.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("--adiabatic").help("-- Enable transform to adiabatic basis.").default_value(false).implicit_value(true);
    program.add_argument("--align").help("-- Align the wavefunction values to the potential.").default_value(false).implicit_value(true);
    program.add_argument("--savewfn").help("-- Save the time evolution of the wavefunction.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    Acorn::QDYN::Options opt = {
        program.get<double>("-f"       ), program.get<double>("-m"), program.get<double>("-p"         ), program.get<double>("-s"     ), program.get<int >("-d"       ),
        program.get<int   >("-i"       ), program.get<int   >("-o"), program.get<bool  >("--adiabatic"), program.get<bool  >("--align"), program.is_used("--optimize" ),
        program.get<bool  >("--savewfn"),
    };

    // return the options and timers
    return {opt, timers};
}

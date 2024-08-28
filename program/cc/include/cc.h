#pragma once

#include "tensor.h"
#include <argparse.hpp>

namespace Acorn {
    namespace CC {
        // options struct
        struct Options {
            std::vector<int> exc, pts; double thresh; int iters; bool linear;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);

        // custom LCCD functions
        namespace LCCD {
            torch::Tensor amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Emsd, const torch::Tensor& T2, int nos);
            double energy(const torch::Tensor& Jmsa, const torch::Tensor& T2, int nos);
        }

        // custom CCD functions
        namespace CCD {
            torch::Tensor amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Emsd, const torch::Tensor& T2, int nos);
            double energy(const torch::Tensor& Jmsa, const torch::Tensor& T2, int nos);
        }

        // custom CCSD functions
        namespace CCSD {
            std::tuple<torch::Tensor, torch::Tensor> amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Fms, const torch::Tensor& Emss, const torch::Tensor& Emsd, const torch::Tensor& T1, const torch::Tensor& T2, int nos);
            double perturbationTriple(const torch::Tensor& Jmsa, const torch::Tensor& Emst, const torch::Tensor& T1, const torch::Tensor& T2, int nos);
            double energy(const torch::Tensor& Jmsa, const torch::Tensor& Fms, const torch::Tensor& T1, const torch::Tensor& T2, int nos);
        }
    }
}

inline std::tuple<Acorn::CC::Options, std::vector<timepoint>> Acorn::CC::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Coupled Clusters Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Number of SCF iterations.").default_value(100).scan<'i', int>();
    program.add_argument("-e", "--excitations").help("-- Excitations to include in the CC calculation.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<int>{1, 2}).scan<'i', int>();
    program.add_argument("-p", "--perturbations").help("-- Perturbations to include in the CC calculation.").nargs(argparse::nargs_pattern::any).default_value(std::vector<int>{}).scan<'i', int>();
    program.add_argument("-t", "--threshold").help("-- Convergence threshold.").default_value(1e-8).scan<'g', double>();
    program.add_argument("--linear").help("-- Use linearized coupled clusters").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    Acorn::CC::Options opt = {program.get<std::vector<int>>("-e"), program.get<std::vector<int>>("-p"), program.get<double>("-t"), program.get<int>("-i"), program.get<bool>("--linear")};

    // return the options and timers
    return {opt, timers};
}

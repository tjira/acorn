#pragma once

#include "tensor.h"
#include <argparse.hpp>
#include <exprtk.hpp>
#include <Eigen/Core>

namespace Acorn {
    namespace Expression {
        // options struct
        struct Options {
            int dim, points; std::string output; std::vector<double> limits; std::vector<std::string> exprs;
        };

        // expression class
        class Expression {
        public:
            Expression(const std::string& exprstr, const std::vector<std::string>& varstr); double eval(const Eigen::VectorXd& r);

        private:
            exprtk::expression<double> expression; Eigen::VectorXd vars;
        };

        // initialization function
        std::tuple<Options, std::vector<timepoint>> initialize(int argc, char** argv);

        // run function
        void run(const Options& opt, std::vector<timepoint>& timers);
    }
}

inline std::tuple<Acorn::Expression::Options, std::vector<timepoint>> Acorn::Expression::initialize(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Expression Evaluator", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<timepoint> timers(2, std::chrono::high_resolution_clock().now());

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-d", "--dimension").help("-- Dimension of the provided expressions.").default_value(1).scan<'i', int>();
    program.add_argument("-e", "--expression").help("-- The expression to evaluate.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<std::string>{"0.5*x^2"});
    program.add_argument("-g", "--grid").help("-- Limits of the evaluation grid.").nargs(2).default_value(std::vector<double>{-16.0, 16.0}).scan<'g', double>();
    program.add_argument("-p", "--points").help("-- Number of points on the grid in each dimension.").default_value(1024).scan<'i', int>();
    program.add_argument("-o", "--output").help("-- Output file path.");

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // initialize the options
    Options opt = {
        program.get<int>("-d"), program.get<int>("-p"), program.get<std::string>("-o"), program.get<std::vector<double>>("-g"), program.get<std::vector<std::string>>("-e")
    };

    // return the options and timers
    return {opt, timers};
}

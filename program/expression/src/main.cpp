#include "expression.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Expression Evaluator", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-e", "--expression").help("-- The expression to evaluate.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<std::string>{"0.5*x^2"});
    program.add_argument("-g", "--grid").help("-- Limits of the evaluation grid.").nargs(2).default_value(std::vector<double>{-16.0, 16.0}).scan<'g', double>();
    program.add_argument("-p", "--points").help("-- Number of points on the grid in each dimension.").default_value(1024).scan<'i', int>();
    program.add_argument("-o", "--output").help("-- Output file path.");

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // print the expression timer label
    std::cout << "EVALUATING THE EXPRESSION: " << std::flush;

    // define the potential matrix
    EigenMatrix<> U(program.get<int>("-p"), program.get<std::vector<std::string>>("-e").size() + 1);

    // calculate the grid spacing
    double dr = (program.get<std::vector<double>>("-g").at(1) - program.get<std::vector<double>>("-g").at(0)) / (program.get<int>("-p") - 1);

    // fill the independent variable column
    for (int i = 0; i < program.get<int>("-p"); i++) U(i, 0) = program.get<std::vector<double>>("-g").at(0) + i * dr;

    // fill the function value matrix columns
    for (int i = 0; i < U.cols() - 1; i++) {
        Expression expr(program.get<std::vector<std::string>>("-e").at(i), {"x"}); U.col(i + 1) = expr.eval(U.col(0));
    }

    // print the elapsed time
    std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // write the expression to disk
    MEASURE("WRITING THE MATRIX:        ", Eigen::Write(program.get<std::string>("-o"), U));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

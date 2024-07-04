#include "expression.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Expression Evaluator", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

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
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the variables
    int dim = program.get<int>("-d"); int points = program.get<int>("-p"); std::vector<double> limits = program.get<std::vector<double>>("-g"); std::vector<std::string> exprs = program.get<std::vector<std::string>>("-e");

    // define the variables vector
    std::vector<std::string> vars(dim);

    // fill the variables vector if the dimension is less than 4
    for (int i = 0; i < dim && dim < 4; i++) vars.at(i) = std::vector<std::string>{"x", "y", "z"}.at(i);

    // fill the variables vector if the dimension is greater than 3
    for (int i = 0; i < dim && dim > 3; i++) vars.at(i) = "x" + std::to_string(i + 1);

    // print the expression timer label
    std::cout << "EVALUATING THE EXPRESSION: " << std::flush;

    // define the expression matrix
    Matrix U((int)std::pow(points, dim), exprs.size() + dim);

    // calculate the grid spacing
    double dr = (limits.at(1) - limits.at(0)) / (points - 1);

    // fill the leftmost independent variable column
    for (int i = 0; i < U.rows(); i++) {
        U.rightCols(exprs.size() + 1).col(0)(i) = limits.at(0) + (i % points) * dr;
    }

    // fill the rest of the independent variable columns
    for (size_t i = 0; i < U.cols() - exprs.size() - 1; i++) {
        for (int j = 0; j < U.rows(); j++) {
            U(j, i) = U(j / (int)std::round(std::pow(points, U.cols() - exprs.size() - i - 1)), U.cols() - exprs.size() - 1);
        }
    }

    // fill the function value matrix columns
    for (int i = 0; i < exprs.size(); i++) {
        Expression expr(exprs.at(i), vars); for (int j = 0; j < U.rows(); j++) U(j, i + dim) = expr.eval(Vector(U.row(j).leftCols(dim)));
    }

    // print the elapsed time
    std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // write the expression to disk
    MEASURE("WRITING THE MATRIX:        ", Eigen::Write(program.get<std::string>("-o"), U));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

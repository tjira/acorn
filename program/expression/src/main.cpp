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

    // fill the variables vector
    for (int i = 0; i < dim && dim < 4; i++) {vars.at(i) = std::vector<std::string>{"x", "y", "z"}.at(i);}

    // print the expression timer label
    std::cout << "EVALUATING THE EXPRESSION: " << std::flush;

    // define the expression matrix
    Matrix U((int)std::pow(points, dim), exprs.size() + dim);

    // calculate the grid spacing
    double dr = (limits.at(1) - limits.at(0)) / (points - 1);

    // fill the independent variable column
    if (dim == 1) {
        for (int i = 0; i < points; i++) {
            U(i, 0) = limits.at(0) + i * dr;
        }
    } else if (dim == 2) {
        for (int i = 0; i < points; i++) {
            for (int j = 0; j < points; j++) {
                U(i * points + j, 0) = limits.at(0) + i * dr;
                U(i * points + j, 1) = limits.at(0) + j * dr;
            }
        }
    } else if (dim == 3) {
        for (int i = 0; i < points; i++) {
            for (int j = 0; j < points; j++) {
                for (int k = 0; k < points; k++) {
                    U(i * points * points + j * points + k, 0) = limits.at(0) + i * dr;
                    U(i * points * points + j * points + k, 1) = limits.at(0) + j * dr;
                    U(i * points * points + j * points + k, 2) = limits.at(0) + k * dr;
                }
            }
        }
    } else {
        throw std::runtime_error(std::to_string(dim) + "-DIMENSIONAL EXPRESSIONS NOT SUPPORTED");
    }

    // fill the function value matrix columns
    for (int i = 0; i < U.cols() - dim; i++) {
        Expression expr(exprs.at(i), vars); for (int j = 0; j < U.rows(); j++) U(j, i + dim) = expr.eval(Vector(U.row(j).leftCols(dim)));
    }

    // print the elapsed time
    std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // write the expression to disk
    MEASURE("WRITING THE MATRIX:        ", Eigen::Write(program.get<std::string>("-o"), U));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

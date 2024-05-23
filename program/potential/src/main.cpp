#include "expression.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-e", "--expression").help("-- The potential energy expression to evaluate.").default_value(std::string("0.5*x^2"));
    program.add_argument("-g", "--grid").help("-- Limits of the evaluation grid.").nargs(2).default_value(std::vector<double>{-16.0, 16.0});
    program.add_argument("-p", "--points").help("-- Number of points on the grid in each dimension.").default_value(1024);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // define the expression and the resulting matrix
    Expression expr(program.get<std::string>("-e"), {"x"}); EigenMatrix<> U(program.get<int>("-p"), 2);

    // calculate the grid spacing
    double dr = (program.get<std::vector<double>>("-g").at(1) - program.get<std::vector<double>>("-g").at(0)) / (program.get<int>("-p") - 1);

    // evaluate the expression on the grid
    for (int i = 0; i < program.get<int>("-p"); i++) {
        U(i, 0) = program.get<std::vector<double>>("-g").at(0) + i * dr; U(i, 1) = expr.eval(U(i, 0));
    }

    // write the potential to disk
    Eigen::Write("U.mat", U);
}

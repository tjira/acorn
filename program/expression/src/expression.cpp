#include "expression.h"

void writeMatrix(const std::string& path, const Eigen::MatrixXd& A) {
    // open the output file and write the dimensions to the header
    std::ofstream file(path); file << A.rows() << " " << A.cols() << "\n" << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < A.rows(); i++, file << "\n") for (int j = 0; j < A.cols(); j++) {
        file << std::setw(20) << A(i, j) << (j < A.cols() - 1 ? " " : "");
    }
}

Acorn::Expression::Expression::Expression(const std::string& exprstr, const std::vector<std::string>& varstr) : vars(varstr.size()) {
    // define the symbol table
    exprtk::symbol_table<double> symbols;

    // add variables and constants to the expression
    for (size_t i = 0; i < varstr.size(); i++) {
        symbols.add_variable(varstr.at(i), vars(i));
    } symbols.add_constants();

    // register the symbol table
    expression.register_symbol_table(symbols);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);
}

double Acorn::Expression::Expression::eval(const Eigen::VectorXd& r) {
    vars = r; return expression.value();
}

void Acorn::Expression::run(const Options& opt, std::vector<timepoint>& timers) {
    // define the variables vector
    std::vector<std::string> vars(opt.dim);

    // fill the variables vector if the dimension is less than 4
    for (int i = 0; i < opt.dim && opt.dim < 4; i++) vars.at(i) = std::vector<std::string>{"x", "y", "z"}.at(i);

    // fill the variables vector if the dimension is greater than 3
    for (int i = 0; i < opt.dim && opt.dim > 3; i++) vars.at(i) = "x" + std::to_string(i + 1);

    // print the expression timer label
    timers.at(1) = std::chrono::high_resolution_clock().now(); std::cout << "EVALUATING THE EXPRESSION: " << std::flush;

    // define the expression matrix
    Eigen::MatrixXd U((int)std::pow(opt.points, opt.dim), opt.exprs.size() + opt.dim);

    // calculate the grid spacing
    double dr = (opt.limits.at(1) - opt.limits.at(0)) / (opt.points - 1);

    // fill the leftmost independent variable column
    for (int i = 0; i < U.rows(); i++) {
        U.rightCols(opt.exprs.size() + 1).col(0)(i) = opt.limits.at(0) + (i % opt.points) * dr;
    }

    // fill the rest of the independent variable columns
    for (size_t i = 0; i < U.cols() - opt.exprs.size() - 1; i++) {
        for (int j = 0; j < U.rows(); j++) {
            U(j, i) = U(j / (int)std::round(std::pow(opt.points, U.cols() - opt.exprs.size() - i - 1)), U.cols() - opt.exprs.size() - 1);
        }
    }

    // fill the function value matrix columns
    for (int i = 0; i < opt.exprs.size(); i++) {
        Acorn::Expression::Expression expr(opt.exprs.at(i), vars); for (int j = 0; j < U.rows(); j++) U(j, i + opt.dim) = expr.eval(Eigen::VectorXd(U.row(j).leftCols(opt.dim)));
    }

    // print the elapsed time
    std::cout << eltime(timers.at(1)) << std::endl;

    // start the timer for writing the matrix
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of the matrix writing
    std::cout << "WRITING THE MATRIX: " << std::flush;

    // write the expression to disk
    writeMatrix(opt.output, U);

    // print the elapsed time for writing the matrix
    std::cout << eltime(timers.at(1)) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << eltime(timers.at(0)) << std::endl;
}

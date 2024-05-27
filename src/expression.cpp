#include "expression.h"

Expression::Expression(const std::string& exprstr, const std::vector<std::string>& varstr) : vars(varstr.size()) {
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

EigenVector<> Expression::eval(const EigenMatrix<>& r) {
    // initialize the result vector
    EigenVector<> result(r.rows());

    // evaluate the expression for each row
    for (int i = 0; i < r.rows(); i++) {
        for (int j = 0; j < r.cols(); j++) {
            vars(j) = r(i, j);
        } result(i) = expression.value();
    }

    // return the result vector
    return result;
}

double Expression::eval(double r) {
    vars(0) = r; return expression.value();
}

#pragma once

#include <Eigen/Core>
#include <exprtk.hpp>

class Expression {
public:
    Expression(const std::string& exprstr, const std::vector<std::string>& varstr); double eval(const Eigen::VectorXd& r);

private:
    exprtk::expression<double> expression; Eigen::VectorXd vars;
};

inline Expression::Expression(const std::string& exprstr, const std::vector<std::string>& varstr) : vars(varstr.size()) {
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

inline double Expression::eval(const Eigen::VectorXd& r) {
    vars = r; return expression.value();
}

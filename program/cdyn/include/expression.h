#pragma once

#include <exprtk.hpp>

class Expression {
public:
    Expression(const std::string& exprstr, const std::vector<std::string>& varstr); double eval(double r);

private:
    exprtk::expression<double> expression; std::vector<double> vars;
};

inline Expression::Expression(const std::string& exprstr, const std::vector<std::string>& varstr) : vars(varstr.size()) {
    // define the symbol table
    exprtk::symbol_table<double> symbols;

    // add variables and constants to the expression
    for (size_t i = 0; i < varstr.size(); i++) {
        symbols.add_variable(varstr.at(i), vars.at(i));
    } symbols.add_constants();

    // register the symbol table
    expression.register_symbol_table(symbols);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);
}

inline double Expression::eval(double r) {
    vars.at(0) = r; return expression.value();
}

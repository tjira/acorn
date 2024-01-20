#include "expression.h"

// include the exprtk
#include <exprtk.hpp>

Vector<> Expression::eval(const Vector<>& r) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x1_exprtk; symbol_table.add_variable("x", x1_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // define the output vector
    Vector<> out(r.size());

    // fill the output vector
    for (int i = 0; i < r.size(); i++) {
        x1_exprtk = r(i); out(i) = expression.value();
    }

    // return the output
    return out;
}

#include "expression.h"

// include the exprtk
#include <exprtk.hpp>

double Expression::eval(double x) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk; symbol_table.add_variable("x", x_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // return the value
    x_exprtk = x; return expression.value();
}

double Expression::eval(double x, double y) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk, y_exprtk; symbol_table.add_variable("x", x_exprtk), symbol_table.add_variable("y", y_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // return the value
    x_exprtk = x, y_exprtk = y; return expression.value();
}

Vector<> Expression::eval(const Vector<>& x) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk; symbol_table.add_variable("x", x_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // define the output vector
    Vector<> out(x.size());

    // fill the output vector
    for (int i = 0; i < x.size(); i++) {
        x_exprtk = x(i); out(i) = expression.value();
    }

    // return the output
    return out;
}

Matrix<> Expression::eval(const Matrix<>& x, const Matrix<>& y) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk, y_exprtk; symbol_table.add_variable("x", x_exprtk), symbol_table.add_variable("y", y_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // define the output vector
    Matrix<> out(x.rows(), x.cols());

    // fill the output vector
    for (int i = 0; i < x.rows(); i++) {
        for (int j = 0; j < x.rows(); j++) {
            x_exprtk = x(i, j), y_exprtk = y(i, j); out(i, j) = expression.value();
        }
    }

    // return the output
    return out;
}

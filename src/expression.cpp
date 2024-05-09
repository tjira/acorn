#include "expression.h"

// include exprtk
#include <exprtk.hpp>

Vector<> Expression::eval(const Vector<>& x) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk; symbol_table.add_variable(variables.at(0), x_exprtk); symbol_table.add_constants();

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
    double x_exprtk, y_exprtk; symbol_table.add_variable(variables.at(0), x_exprtk), symbol_table.add_variable(variables.at(1), y_exprtk); symbol_table.add_constants();

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

double Expression::get(double x) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk; symbol_table.add_variable(variables.at(0), x_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // return the value
    x_exprtk = x; return expression.value();
}

double Expression::get(double x, double y) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk, y_exprtk; symbol_table.add_variable(variables.at(0), x_exprtk), symbol_table.add_variable(variables.at(1), y_exprtk); symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // return the value
    x_exprtk = x, y_exprtk = y; return expression.value();
}

double Expression::get(const Vector<>& r) const {
    // define the exprtk symbol table and expression
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;

    // add variables and constants to the expression
    double x_exprtk, y_exprtk;
    if (r.size() > 0) symbol_table.add_variable(variables.at(0), x_exprtk);
    if (r.size() > 1) symbol_table.add_variable(variables.at(1), y_exprtk);
    symbol_table.add_constants();

    // register the symbol table with the expression
    expression.register_symbol_table(symbol_table);

    // parse the expression
    exprtk::parser<double>().compile(exprstr, expression);

    // return the value
    x_exprtk = r.size() > 0 ? r(0) : 0, y_exprtk = r.size() > 1 ? r(1) : 0; return expression.value();
}

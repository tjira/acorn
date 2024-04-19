#pragma once

#include "eigen.h"
#include "numpy.h"

class Expression {
public:
    // constructor
    Expression(const std::string& expr, const std::vector<std::string>& variables) : variables(variables), exprstr(expr) {}

    // methods
    double get(double x) const; double get(double x, double y) const; double get(const Vector<>& r) const;
    Vector<> eval(const Vector<>& x) const; Matrix<> eval(const Matrix<>& x, const Matrix<>& y) const;

private:
    std::vector<std::string> variables; std::string exprstr;
};

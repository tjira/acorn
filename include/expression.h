#pragma once

#include "eigen.h"
#include "numpy.h"

class Expression {
public:
    // constructor
    Expression(const std::string& expr) : exprstr(expr) {}

    // methods
    double eval(double x) const; double eval(double x, double y) const; Vector<> eval(const Vector<>& x) const; Matrix<> eval(const Matrix<>& x, const Matrix<>& y) const;

private:
    std::string exprstr;
};

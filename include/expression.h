#pragma once

#include "eigen.h"

class Expression {
public:
    // constructor
    Expression(const std::string& expr) : exprstr(expr) {}

    // methods
    double eval(double r) const; Vector<> eval(const Vector<>& r) const;

private:
    std::string exprstr;
};

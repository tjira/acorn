#pragma once

#include "linalg.h"
#include <exprtk.hpp>

class Expression {
public:
    Expression(const std::string& exprstr, const std::vector<std::string>& varstr);

    // evaluation functions
    double eval(double r); double eval(const Vector& r);

private:
    exprtk::expression<double> expression; Vector vars;
};

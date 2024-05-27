#pragma once

#include "linalg.h"
#include <exprtk.hpp>

class Expression {
public:
    Expression(const std::string& exprstr, const std::vector<std::string>& varstr);

    // evaluation functions
    EigenVector<> eval(const EigenMatrix<>& r); double eval(double r);

private:
    exprtk::expression<double> expression; EigenVector<> vars;
};

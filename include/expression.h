#pragma once

#include "exprtk.h"
#include "ptable.h"
#include "system.h"

class Expression {
public:
    Expression(const std::string& exprstr, const std::vector<std::string>& varstr);

    // evaluation functions
    double eval(double r);

private:
    exprtk::expression<double> expression; EigenVector<> vars;
};

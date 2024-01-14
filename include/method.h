#pragma once

#include "integral.h"

class Method {
public:
    struct Result {
        struct RestrictedHartreeFock {
            double E; Vector<> eps; Matrix<> C, D;
        } rhf;
    };
public:
    virtual Result run(const System&, const Integrals&, bool) const = 0;
};

#include <nlohmann/json.hpp>

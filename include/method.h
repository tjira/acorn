#pragma once

#include "integral.h"

class Method {
public:
    struct Result {
        struct RestrictedHartreeFock {
            double E; Vector<> eps; Matrix<> C, D;
        } rhf;
        struct RestrictedMollerPlesset {
            double Ecorr = 0;
        } rmp;
        double Etot; Matrix<> G;
    };
public:
    virtual Result gradient(const System&, const Integrals&, Result, bool print = true) const;
    virtual Result run(const System&, const Integrals&, Result, bool print = true) const = 0;
};

#include <nlohmann/json.hpp>

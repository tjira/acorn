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
    // general methods
    Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const;

    // virtual functions
    virtual double energy(const System&) const = 0;
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
};

#include <nlohmann/json.hpp>

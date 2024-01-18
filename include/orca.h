#pragma once

#include "method.h"

class Orca : public Method<Orca> {
    friend class Method<Orca>;
public:
    struct Options {Options(){};
        // structures
        struct Dynamics {int iters=100; double step=1.0; std::string output="trajectory.xyz";} dynamics={};
        struct Gradient {double step=1e-5;} gradient={}; struct Hessian {double step=1e-5;} hessian={};

        // variables
        std::filesystem::path folder; std::string interface;
    };
public:
    // constructors and destructors
    Orca(const Options& opt = {}) : opt(opt) {}

    // overriden virtual derivatives
    Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const override;

    // overriden virtual methods
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;
    Result run(const System& system, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

#pragma once

#include "restrictedmethod.h"

class RestrictedHartreeFock : public RestrictedMethod<RestrictedHartreeFock> {
    friend class Method<RestrictedHartreeFock>;
public:
    struct Options {
        // structures
        struct Gradient {double step;} gradient; struct Hessian {double step;} hessian;
        struct Dynamics {int iters; double step; std::string output;} dynamics;

        // variables
        int maxiter; double thresh;
    };
public:
    // constructors and destructors
    RestrictedHartreeFock(const Options& opt) : opt(opt) {}

    // overriden virtual methods
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;
    Result run(const System& system, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

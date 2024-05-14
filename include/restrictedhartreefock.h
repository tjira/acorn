#pragma once

#include "restrictedmethod.h"

class RestrictedHartreeFock : public RestrictedMethod<RestrictedHartreeFock> {
    friend class Method<RestrictedHartreeFock>;
public:
    struct Options {Options(){};
        // gradient and hessian structrures
        struct Gradient {double step;} gradient={}; struct Hessian {double step;} hessian={};

        // dynamics structure
        struct Dynamics {
            int iters; double step;
        } dynamics={};

        // variables
        int maxiter; double thresh;
    };
public:
    // constructors and destructors
    RestrictedHartreeFock(const Options& opt = {}) : opt(opt) {}

    // overriden virtual methods
    Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const override;
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;

private:
    std::tuple<Result, Integrals> run(const System& system, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

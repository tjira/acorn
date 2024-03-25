#pragma once

#include "restrictedmethod.h"

class RestrictedHartreeFock : public RestrictedMethod<RestrictedHartreeFock> {
    friend class Method<RestrictedHartreeFock>;
public:
    struct Options {Options(){};
        // gradient and hessian structrures
        struct Gradient {double step=1e-5;} gradient={}; struct Hessian {double step=1e-5;} hessian={};

        // dynamics structure
        struct Dynamics {
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0;
        } dynamics={};

        // variables
        int maxiter=100; double thresh=1e-8;
    };
public:
    // constructors and destructors
    RestrictedHartreeFock(const Options& opt = {}) : opt(opt) {}

    // overriden virtual methods
    Result gradient(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;

private:
    std::tuple<Result, Integrals> run(const System& system, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

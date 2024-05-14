#pragma once

#include "unrestrictedmethod.h"

class UnrestrictedHartreeFock : public UnrestrictedMethod<UnrestrictedHartreeFock> {
    friend class Method<UnrestrictedHartreeFock>;
public:
    struct Options {Options(){};
        // gradient and hessian structrures
        struct Gradient {double step;} gradient={}; struct Hessian {double step;} hessian={};

        // dynamics structure
        struct Dynamics {
            int iters; double step; std::string folder;
        } dynamics={};

        // variables
        int maxiter; double thresh;
    };
public:
    // constructors and destructors
    UnrestrictedHartreeFock(const Options& opt = {}) : opt(opt) {}

    // overriden virtual methods
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;

private:
    std::tuple<Result, Integrals> run(const System& system, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

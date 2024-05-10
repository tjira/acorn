#pragma once

#include "restrictedmollerplesset.h"
#include "determinant.h"

class RestrictedConfigurationInteraction : public RestrictedMethod<RestrictedConfigurationInteraction> {
    friend class Method<RestrictedConfigurationInteraction>;
public:
    struct Options {Options(){};
        // gradient and hessian structures
        struct Gradient {double step=1e-5;} gradient; struct Hessian {double step=1e-5;} hessian;

        // dynamics structure
        struct Dynamics {
            int iters=100; double step=1.0;
        } dynamics={};

        // RCI options
        std::vector<int> excitations; int nstate=2, state=1;
    };
public:
    // constructors and destructors
    RestrictedConfigurationInteraction(const RestrictedHartreeFock::Options& rhfopt = {}, const Options& opt = {}) : opt(opt), rhfopt(rhfopt) {}

    // overriden virtual methods
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    std::tuple<Result, Integrals> run(const System& system, Result res, bool print = true) const override;

private:
    Options opt; RestrictedHartreeFock::Options rhfopt;
};

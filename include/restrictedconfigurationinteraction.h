#pragma once

#include "restrictedmollerplesset.h"
#include "determinant.h"

class RestrictedConfigurationInteraction : public RestrictedMethod<RestrictedConfigurationInteraction> {
    friend class Method<RestrictedConfigurationInteraction>;
public:
    struct Options {Options(){};
        // structure of the options
        struct Gradient {double step=1e-5;} gradient; struct Hessian {double step=1e-5;} hessian;
        struct Dynamics {int iters=100; double step=1.0; std::string folder;} dynamics;
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

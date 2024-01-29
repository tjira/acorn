#pragma once

#include "restrictedmollerplesset.h"

class RestrictedConfigurationInteraction : public RestrictedMethod<RestrictedConfigurationInteraction> {
    friend class Method<RestrictedConfigurationInteraction>;
public:
    struct Options {Options(){};
        // structure of the options
        struct Dynamics {int iters=100; double step=1.0; std::string output="trajectory.xyz";} dynamics;
        struct Gradient {double step=1e-5;} gradient; struct Hessian {double step=1e-5;} hessian;
    };
public:
    // constructors and destructors
    RestrictedConfigurationInteraction(const RestrictedHartreeFock::Options& rhfopt, const Options& opt) : opt(opt), rhfopt(rhfopt) {}

    // overriden virtual methods
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    Result run(const System& system, Result res, bool print = true) const override;

private:
    Options opt; RestrictedHartreeFock::Options rhfopt;
};

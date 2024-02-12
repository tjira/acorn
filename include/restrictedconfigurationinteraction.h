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
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0; std::string folder;
        } dynamics={};

        // RCI options
        std::vector<int> excitations;
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

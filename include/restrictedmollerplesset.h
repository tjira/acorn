#pragma once

#include "restrictedhartreefock.h"
#include "transform.h"

class RestrictedMollerPlesset : public RestrictedMethod<RestrictedMollerPlesset> {
    friend class Method<RestrictedMollerPlesset>;
public:
    struct Options {Options(){};
        // gradient and hessian structures
        struct Gradient {double step=1e-5;} gradient; struct Hessian {double step=1e-5;} hessian;

        // dynamics structure
        struct Dynamics {
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0;
        } dynamics={};

        // values of simple options
        int order=2;
    };
public:
    // constructors and destructors
    RestrictedMollerPlesset(const RestrictedHartreeFock::Options& rhfopt = {}, const Options& opt = {}) : opt(opt), rhfopt(rhfopt) {}

    // overriden virtual methods
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    std::tuple<Result, Integrals> run(const System& system, Result res, bool print = true) const override;

private:
    Options opt; RestrictedHartreeFock::Options rhfopt;
};

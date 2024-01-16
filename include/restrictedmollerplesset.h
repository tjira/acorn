#pragma once

#include "restrictedhartreefock.h"
#include "transform.h"

class RestrictedMollerPlesset : public RestrictedMethod<RestrictedMollerPlesset> {
    friend class Method<RestrictedMollerPlesset>;
public:
    struct Options {
        // structure of the options
        struct Gradient {double step;} gradient; struct Hessian {double step;} hessian;
        struct Dynamics {int iters; double step; std::string output;} dynamics;

        // values of simple options
        int order;
    };
public:
    // constructors and destructors
    RestrictedMollerPlesset(const Options& opt, const RestrictedHartreeFock::Options& rhfopt) : opt(opt), rhfopt(rhfopt) {}

    // methods
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;
    Result run(const System& system, Result res, bool print = true) const override;

private:
    Options opt; RestrictedHartreeFock::Options rhfopt;
};

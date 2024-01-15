#pragma once

#include "restrictedhartreefock.h"
#include "transform.h"

class RestrictedMollerPlesset : public RestrictedMethod<RestrictedMollerPlesset> {
    friend class Method<RestrictedMollerPlesset>;
public:
    struct Options {
        // parent method options
        RestrictedHartreeFock::Options rhfopt;

        // structure of the options
        struct Gradient {double step;} gradient; struct Hessian {double step;} hessian;
        struct Dynamics {int iters; double step; std::string output;} dynamics;

        // values of simple options
        int order = rmpopt.at("order");
    };
public:
    // constructors and destructors
    RestrictedMollerPlesset(const Options& opt) : opt(opt) {}

    // methods
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;
    Result run(const System& system, Result res, bool print = true) const override;

private:
    Options opt;
};

// option structures loaders
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Hessian, step);

// options loader
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, dynamics, gradient, hessian, rhfopt, order);

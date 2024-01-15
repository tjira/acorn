#pragma once

#include "restrictedhartreefock.h"
#include "transform.h"

class RestrictedMollerPlesset : public RestrictedMethod {
public:
    struct Options {
        // structure of the options
        RestrictedHartreeFock::Options rhfopt;
        Method::Options::Gradient gradient;
        Method::Options::Hessian hessian;

        // values of simple options
        int order = rmpopt.at("order");
    };
public:
    // constructors and destructors
    RestrictedMollerPlesset(const Options& opt) : RestrictedMethod({opt.gradient, opt.hessian}), opt(opt) {}

    // methods
    double energy(const System& system, Result res) const override;
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    Options opt;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, gradient, hessian, rhfopt, order);

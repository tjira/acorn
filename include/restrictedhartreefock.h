#pragma once

#include "restrictedmethod.h"

class RestrictedHartreeFock : public RestrictedMethod {
public:
    struct Options {
        // structure of the options
        Method::Options::Gradient gradient;
        Method::Options::Hessian hessian;

        // values of simple options
        double thresh = rhfopt.at("thresh");
        int maxiter = rhfopt.at("maxiter");
    };
public:
    // constructors and destructors
    RestrictedHartreeFock(const Options& opt) : RestrictedMethod({opt.gradient, opt.hessian}), opt(opt) {}

    // overriden virtual methods
    double energy(const System& system, Result res = {}) const override;
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options, gradient, hessian, maxiter, thresh);

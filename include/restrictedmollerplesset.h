#pragma once

#include "restrictedhartreefock.h"
#include "transform.h"

class RestrictedMollerPlesset : public RestrictedMethod {
public:
    struct Options {
        RestrictedHartreeFock::Options rhfopt = {}; int order = 2;
    };
public:
    // constructors and destructors
    RestrictedMollerPlesset(const Options& opt) : opt(opt) {}

    // methods
    double energy(const System& system) const override;
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    Options opt;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, rhfopt, order);

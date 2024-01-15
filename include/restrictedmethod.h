#pragma once

#include "method.h"

class RestrictedMethod : public Method {
public:
    // constructors and destructors
    RestrictedMethod(const Method::Options& opt) : Method(opt) {}

    // virtual functions
    virtual double energy(const System&, Result) const = 0;
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
};

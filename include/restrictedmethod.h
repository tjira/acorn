#pragma once

#include "method.h"

class RestrictedMethod : public Method {
public:
    virtual double energy(const System&) const = 0;
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
};

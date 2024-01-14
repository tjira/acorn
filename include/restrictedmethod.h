#pragma once

#include "method.h"

class RestrictedMethod : public Method {
public:
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
};

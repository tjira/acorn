#pragma once

#include "method.h"

template <class M>
class RestrictedMethod : public Method<M> {
public:
    // virtual functions
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
    virtual Result run(const System&, Result, bool) const = 0;
};

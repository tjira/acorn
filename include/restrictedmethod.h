#include "method.h"

class RestrictedMethod : public Method {
public:
    virtual Result run(const System&, const Integrals&, bool) const = 0;
};

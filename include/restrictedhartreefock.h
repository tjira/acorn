#include "restrictedmethod.h"

class RestrictedHartreeFock : public RestrictedMethod {
public:
    struct Options {
        Options() = default; int maxiter = 100; double thresh = 1e-8;
    };
public:
    // constructors and destructors
    RestrictedHartreeFock(const Options& opt) : opt(opt) {}

    // methods
    Result run(const System& system, const Integrals& ints, Result res = {}, bool print = true) const override;

private:
    Options opt;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options, maxiter, thresh);

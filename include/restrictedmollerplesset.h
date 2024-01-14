#include "restrictedmethod.h"
#include "transform.h"

class RestrictedMollerPlesset : public RestrictedMethod {
public:
    struct Options {
        int order = 2;
    };
public:
    // constructors and destructors
    RestrictedMollerPlesset(const Options& opt) : opt(opt) {}

    // methods
    Result run(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    Options opt;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, order);

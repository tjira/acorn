#include "restrictedhartreefock.h"
#include "population.h"
#include "printer.h"

static std::function<Vector<>()> result = []() {Vector<> res(3);
    res << -0.33053075785438, 0.16526569728403, 0.16526506057035;
return res;}; constexpr double precision = 1e-6;

int test_rhf_mulliken(int, char**) {
    // create the molecule stream
    std::ifstream mstream("../example/molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "STO-3G"); Integrals ints(true);

    // define the HF options
    RestrictedHartreeFock::Options rhfopt; rhfopt.maxiter = 100; rhfopt.thresh = 1e-8;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = RestrictedHartreeFock(rhfopt).run(system, ints, {}, false);

    // run the mullikem population analysis
    Vector<> mulliken = Population::Mulliken(system, ints, res.rhf.D);

    // set the printer flags
    std::cout << std::fixed << std::setprecision(14);

    // print the population and the differences from the reference result
    Printer::Print(mulliken, "MULLIKEN POPULATIONS"); std::cout << std::endl; Printer::Print((mulliken - result()).array().abs(), "MULLIKEN POPULATION DIFFERENCES");

    // return the status of the test
    return ((mulliken - result()).array().abs() < precision).all() ? 0 : 1;
}

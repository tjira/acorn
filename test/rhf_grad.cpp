#include "restrictedhartreefock.h"
#include "population.h"
#include "printer.h"

static std::function<Matrix<>()> result = []() {Matrix<> res(3, 3);
    res <<  0.00000206715938, -0.00000084810144,  0.00000127834436,
           -0.00000144284456, -0.00000072895739, -0.00000051094623,
           -0.00000062431486,  0.00000157705885, -0.00000076739812;
return res;}; constexpr double precision = 1e-6;

int test_rhf_grad(int, char**) {
    // create the molecule stream
    std::ifstream mstream("../example/molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "STO-3G"); Integrals ints(true);

    // define the HF options
    RestrictedHartreeFock::Options rhfopt; rhfopt.maxiter = 100; rhfopt.thresh = 1e-8;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // calculate atomic integrals derivatives
    ints.dS = Integral::dOverlap(system); ints.dT = Integral::dKinetic(system);
    ints.dV = Integral::dNuclear(system); ints.dJ = Integral::dCoulomb(system);

    // create the Hartree-Fock object
    RestrictedHartreeFock rhf(rhfopt);

    // run the restricted Hartree-Fock calculation
    Result res = rhf.run(system, ints, {}, false);

    // calculate the analytical gradient
    res = rhf.gradient(system, ints, res);

    // set the printer flags
    std::cout << std::fixed << std::setprecision(14);

    // print the population and the differences from the reference result
    Printer::Print(res.rhf.G, "GRADIENT"); std::cout << std::endl; Printer::Print((res.rhf.G - result()).array().abs(), "GRADIENT DIFFERENCES");

    // return the status of the test
    return ((res.rhf.G - result()).array().abs() < precision).all() ? 0 : 1;
}

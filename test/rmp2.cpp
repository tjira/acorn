#include "restrictedmollerplesset.h"

constexpr double result = -75.00485494762017, precision = 1e-6;

int test_rmp2(int, char**) {
    // create the molecule stream
    std::ifstream mstream("../example/molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "STO-3G"); Integrals ints(true);

    // define the HF and MP2 options
    RestrictedHartreeFock::Options rhfopt; rhfopt.maxiter = 100; rhfopt.thresh = 1e-8;
    RestrictedMollerPlesset::Options rmpopt; rmpopt.order = 2;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = RestrictedHartreeFock(rhfopt).run(system, ints, {}, false);

    // transform the coulomb integrals to the MO basis
    ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C);

    // perform the RMP2 calculation
    res = RestrictedMollerPlesset(rhfopt, rmpopt).run(system, ints, res, false);

    // print the total energy and the difference from the reference result
    std::printf("ENERGY: %.14f, DIFFERENCE: %.3e\n", res.Etot, std::abs(res.Etot - result));

    // return the status of the test
    return std::abs(res.Etot - result) < precision ? 0 : 1;
}

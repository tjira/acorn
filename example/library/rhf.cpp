#include "acorn.h"

int main() {
    // create the molecule stream
    std::ifstream mstream("../molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "sto-3g"); Integrals ints(true);

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = RestrictedHartreeFock().run(system, ints, {}, false);

    // print the energy
    std::cout << std::setprecision(14) << res.Etot << std::endl;
}

#include "acorn.h"

int main() {
    // create the molecule stream
    std::ifstream mstream("../molecule/water.xyz");

    // create the molecule
    System system(mstream, "sto-3g");

    // create the Integrals and the RHF object
    Integrals ints; RestrictedHartreeFock rhf;

    // calculate all the atomic integrals
    Integral::Initialize();
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);
    Integral::Finalize();

    // run the restricted Hartree-Fock calculation
    Result res = rhf.run(system, ints, {}, false);

    // print the energy
    std::cout << std::setprecision(14) << res.Etot << std::endl;
}

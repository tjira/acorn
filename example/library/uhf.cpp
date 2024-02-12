#include "acorn.h"

int main(int argc, char** argv) {
    // get executable path for the executable to run from anywhere
    auto path = std::filesystem::weakly_canonical(argv[0]).parent_path();

    // set the basis folder path to the project root
    setenv("LIBINT_DATA_PATH", path.parent_path().parent_path().c_str(), true);

    // create the molecule stream
    std::ifstream mstream(path / "../molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "sto-3g", 1, 2); Integrals ints(true);

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = UnrestrictedHartreeFock().run(system, ints, {}, false);

    // print the total energy
    std::printf("%.8f\n", res.Etot);
}

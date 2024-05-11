#include "unrestrictedhartreefock.h"

constexpr double result = -74.66228146841445, precision = 1e-6;

int main(int, char** argv) {
    // get executable path for the executable
    auto path = std::filesystem::weakly_canonical(argv[0]).parent_path();

    // create the molecule stream
    std::ifstream mstream(path / "../../example/molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "STO-3G", 1, 2); Integrals ints(true);

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = UnrestrictedHartreeFock().run(system, ints, {}, false);

    // print the total energy and the difference from the reference result
    std::printf("ENERGY: %.14f, DIFFERENCE: %.3e\n", res.Etot, std::abs(res.Etot - result));

    // return the status of the test
    return std::abs(res.Etot - result) < precision ? 0 : 1;
}

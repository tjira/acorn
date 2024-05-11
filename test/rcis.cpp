#include "restrictedconfigurationinteraction.h"

constexpr double result = -74.96590121728482, precision = 1e-6;

int main(int, char** argv) {
    // get executable path for the executable
    auto path = std::filesystem::weakly_canonical(argv[0]).parent_path();

    // create the molecule stream
    std::ifstream mstream(path / "../../example/molecule/water.xyz");

    // create the molecule and integral container
    System system(mstream, "STO-3G"); Integrals ints(true);

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = RestrictedHartreeFock().run(system, ints, {}, false);

    // transform the one-electron integrals to the MS basis
    ints.Tms = Transform::SingleSpin(ints.T, res.rhf.C);
    ints.Vms = Transform::SingleSpin(ints.V, res.rhf.C);

    // transform the coulomb integrals to the MS basis
    ints.Jms = Transform::CoulombSpin(ints.J, res.rhf.C);

    // create the CI options
    RestrictedConfigurationInteraction::Options opt; opt.excitations = {1};

    // perform the CISD calculation
    res = RestrictedConfigurationInteraction(RestrictedHartreeFock::Options(), opt).run(system, ints, res, false);

    // print the total energy and the difference from the reference result
    std::printf("ENERGY: %.14f, DIFFERENCE: %.3e\n", res.Etot, std::abs(res.Etot - result));

    // return the status of the test
    return std::abs(res.Etot - result) < precision ? 0 : 1;
}

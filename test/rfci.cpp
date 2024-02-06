#include "acorn.h"

int main(int argc, char** argv) {
    // throw an error if the number of arguments is incorrect
    if (argc != 3) throw std::runtime_error("INCORRECT NUMBER OF ARGUMENTS");

    // get executable path for the executable to run from anywhere
    auto path = std::filesystem::weakly_canonical(argv[0]).parent_path();

    // create the molecule stream
    std::ifstream mstream(path / "../../example/molecule" / (argv[1] + std::string(".xyz")));

    // create the molecule and integral container
    System system(mstream, argv[2]); Integrals ints(true);

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = RestrictedHartreeFock().run(system, ints, {}, false);

    // transform the one electron integrals to the MS basis
    ints.Tms = Transform::SingleSpin(ints.T, res.rhf.C);
    ints.Vms = Transform::SingleSpin(ints.V, res.rhf.C);

    // transform the coulomb integrals to the MS basis
    ints.Jms = Transform::CoulombSpin(ints.J, res.rhf.C);

    // perform the FCI calculation
    res = RestrictedConfigurationInteraction().run(system, ints, res, false);

    // print the total energy
    std::printf("%.8f\n", res.Etot);
}

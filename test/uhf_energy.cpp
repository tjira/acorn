#include "../bin/acorn.h"

int main(int argc, char** argv) {
    // throw an error if the number of arguments is incorrect
    if (argc != 5) throw std::runtime_error("INCORRECT NUMBER OF ARGUMENTS");

    // get executable path for the executable to run from anywhere
    auto path = std::filesystem::weakly_canonical(argv[0]).parent_path();

    // create the molecule stream
    std::ifstream mstream(path / "../../example/molecule" / (argv[1] + std::string(".xyz")));

    // create the molecule and integral container
    System system(mstream, argv[2], std::stoi(argv[3]), std::stoi(argv[4])); Integrals ints(true);

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system); ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system); ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock calculation
    Result res = UnrestrictedHartreeFock().run(system, ints, {}, false);

    // print the total energy
    std::printf("%.8f\n", res.Etot);
}

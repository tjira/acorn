#include "restrictedhartreefock.h"
#include <argparse/argparse.hpp>

#include "default.h" // default options

int main(int argc, char** argv) {
    // create the argument parser
    argparse::ArgumentParser program("acorn", "0.0.1", argparse::default_arguments::none);

    // add the command line arguments arguments
    program.add_argument("input").help("-- Input file to specify the calculation.");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads.").default_value(1).scan<'i', int>();

    // parse the command line arguments and print help if requested
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set the number of threads and defne results
    nthread = program.get<int>("-n"); Method::Result res;

    // open the input file, parse the input, create integrals container and create the molecule
    std::ifstream istream(program.get("input")); nlohmann::json input = nlohmann::json::parse(istream); istream.close(); Integrals ints;
    std::ifstream mstream(input.at("molecule").at("file")); System system(mstream, input.at("molecule").at("basis")); mstream.close();

    // print program name and compiler version
    std::printf("QUANTUM ACORN (GCC %d.%d.%d, ", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);

    // print library versions
    std::printf("LIBINT %d.%d.%d, ", LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);
    std::printf("EIGEN %d.%d.%d)\n", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);

    // calculate all the atomic integrals
    libint2::initialize(); ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system); libint2::finalize();

    // path the input json file and apply defaults
    if (input.contains("rhf")) rhfopt.merge_patch(input.at("rhf"));

    if (input.contains("rhf")) {
        // initialize the restricted Hartree-Fock method
        RestrictedHartreeFock rhf(rhfopt.template get<RestrictedHartreeFock::Options>());

        // run the restricted Hartree-Fock method
        res = rhf.run(system, ints);

        // print the final SCF energy
        std::printf("\nSCF ENERGY: %.14f\n", res.rhf.E + system.repulsion());
    }
}

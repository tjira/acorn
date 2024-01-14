#include "restrictedmollerplesset.h"
#include "restrictedhartreefock.h"
#include <argparse/argparse.hpp>

// transform algorithms
#include "transform.h"

int main(int argc, char** argv) {
    // set path to the basis function folder if not set
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        setenv("LIBINT_DATA_PATH", std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path().c_str(), true);
    }

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

    // set the number of threads and defne results and integrals
    nthread = program.get<int>("-n"); Method::Result res; Integrals ints;

    // open the input file, parse the input, create integrals container and create the molecule
    std::ifstream istream(program.get("input")); nlohmann::json input = nlohmann::json::parse(istream); istream.close();

    // get the path of the system file
    auto syspath = std::filesystem::path(program.get("input")).parent_path().string() / std::filesystem::path(input.at("molecule").at("file"));

    // create the system from the system file
    std::ifstream mstream(syspath); System system(mstream, input.at("molecule").at("basis")); mstream.close();

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
    if (input.contains("mp")) mpopt.merge_patch(input.at("mp"));

    if (input.contains("rhf")) {
        // initialize the restricted Hartree-Fock method
        RestrictedHartreeFock rhf(rhfopt.template get<RestrictedHartreeFock::Options>());

        // run the restricted Hartree-Fock method
        res = rhf.run(system, ints);

        // print the final SCF energy
        std::printf("\nSCF ENERGY: %.14f\n", res.Etot);

        if (rhfopt.at("gradient")) {
            // calculate the gradient
            res = rhf.gradient(system, ints, res);

            // print the final gradient
            std::cout << "\nGRADIENT:\n" << res.G << "\n\n" << "GRADIENT NORM: " << res.G.norm() << std::endl;
        }

        if (input.contains("mp")) {
            // initialize the restricted Moller-Plesset method
            RestrictedMollerPlesset rmp(mpopt.template get<RestrictedMollerPlesset::Options>());

            // run the restricted Moller-Plesset method
            res = rmp.run(system, ints, res);

            // print the final MP2 energy
            std::printf("\nMP2 ENERGY: %.14f\n", res.Etot);
        }
    }
}

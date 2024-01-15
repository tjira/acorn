#include "restrictedmollerplesset.h"
#include "restrictedhartreefock.h"
#include <argparse/argparse.hpp>

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

    // patch the input json file and apply defaults
    if (input.contains("rhf")) rhfopt.merge_patch(input.at("rhf"));
    if (input.contains("rmp")) rmpopt.merge_patch(input.at("rmp"));

    // pass the HF options to the post-HF methods
    rmpopt["rhfopt"] = rhfopt.template get<RestrictedHartreeFock::Options>();

    if (input.contains("rhf")) {
        // perform the RHF calculation, assign it to the result struct and print the results
        res = RestrictedHartreeFock(rhfopt.template get<RestrictedHartreeFock::Options>()).run(system, ints);
        std::printf("\nRHF ENERGY: %.14f\n", res.Etot);

        if (libint2::initialize(); input.at("rhf").contains("gradient")) {
            // calculate the RHF gradient, assign it to the result struct and print the results
            res = RestrictedHartreeFock(rhfopt.template get<RestrictedHartreeFock::Options>()).gradient(system, ints, res);
            std::cout << "\nRHF GRADIENT:\n" << res.G << "\n" << "RHF GRADIENT NORM: " << res.G.norm() << std::endl;
        } libint2::finalize();

        if (libint2::initialize(); input.at("rhf").contains("hessian")) {
            // calculate the RHF gradient, assign it to the result struct and print the results
            res = RestrictedHartreeFock(rhfopt.template get<RestrictedHartreeFock::Options>()).hessian(system, ints, res);
            std::cout << "\nRHF HESSIAN:\n" << res.H << "\n" << "RHF HESSIAN NORM: " << res.H.norm() << std::endl;
            std::cout << "\nRHF FREQUENCIES:\n" << Method::frequency(system, res.H) << std::endl;
        } libint2::finalize();

        if (input.contains("rmp")) {
            // initialize the restricted RMP method, run it, assign it to the result struct and print the results
            RestrictedMollerPlesset rmp(rmpopt.template get<RestrictedMollerPlesset::Options>());
            ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C); res = rmp.run(system, ints, res);
            std::printf("\nRMP2 ENERGY: %.14f\n", res.Etot);

            if (libint2::initialize(); input.at("rmp").contains("gradient")) {
                // calculate the RMP gradient, assign it to the result struct and print the results
                res = rmp.gradient(system, ints, res); std::cout << "\nRMP2 GRADIENT:\n" << res.G << "\n" << "RMP2 GRADIENT NORM: " << res.G.norm() << std::endl;
            } libint2::finalize();

            if (libint2::initialize(); input.at("rmp").contains("hessian")) {
                // calculate the RHF gradient, assign it to the result struct and print the results
                res = rmp.hessian(system, ints, res); std::cout << "\nRMP HESSIAN:\n" << res.H << "\n" << "RMP HESSIAN NORM: " << res.H.norm() << std::endl;
                std::cout << "\nRMP FREQUENCIES:\n" << Method::frequency(system, res.H) << std::endl;
            } libint2::finalize();
        }
    }
}

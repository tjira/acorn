// include the methods
#include "restrictedmollerplesset.h"
#include "restrictedhartreefock.h"

// model system and solver
#include "modelsolver.h"

// interfaces
#include "orca.h"

// include default options
#include "default.h"

// include parsing and json libraries
#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>

// option structures loaders for RHF
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Hessian, step);

// option structures loaders for RMP
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Hessian, step);

// option structures loaders for ORCA
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Hessian, step);

// option structures loaders for model method
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::Options::Dynamics, iters, step, output);

// option loaders
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options, dynamics, gradient, hessian, maxiter, thresh);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::Options, dynamics, real, step, iters, nstate, thresh, optimize);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, dynamics, gradient, hessian, order);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options, dynamics, interface, folder);

int main(int argc, char** argv) {
    // create the argument parser and the program timer
    argparse::ArgumentParser program("acorn", "0.0.1", argparse::default_arguments::none); auto programtimer = Timer::Now();

    // add the command line arguments arguments
    program.add_argument("input").help("-- Input file to specify the calculation.");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads.").default_value(1).scan<'i', int>();

    // parse the command line arguments and print help if requested
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set path to the basis function folder if not set
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        setenv("LIBINT_DATA_PATH", std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path().c_str(), true);
    } nthread = program.get<int>("-n"); std::cout << std::fixed << std::setprecision(14);

    // get the path of the input
    auto inputpath = std::filesystem::current_path() / std::filesystem::path(program.get("input")).parent_path();


    // print program name and compiler version
    std::printf("QUANTUM ACORN (GCC %d.%d.%d, ", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);

    // print library versions
    std::printf("EIGEN %d.%d.%d, ", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
    std::printf("LIBINT %d.%d.%d)\n", libint2::major(), libint2::minor(), libint2::micro());

    // open the input file and parse the input
    std::ifstream istream(program.get("input")); nlohmann::json input = nlohmann::json::parse(istream); istream.close();

    // parse the default options
    auto rmpopt = nlohmann::json::parse(rmpoptstr), orcaopt = nlohmann::json::parse(orcaoptstr);
    auto intopt = nlohmann::json::parse(intoptstr), rhfopt = nlohmann::json::parse(rhfoptstr);
    auto msvopt = nlohmann::json::parse(msvoptstr);

    // patch the input json file and apply defaults
    if (input.contains("integral")) intopt.merge_patch(input.at("integral"));
    if (input.contains("solve")) msvopt.merge_patch(input.at("solve"));
    if (input.contains("orca")) orcaopt.merge_patch(input.at("orca"));
    if (input.contains("rhf")) rhfopt.merge_patch(input.at("rhf"));
    if (input.contains("rmp")) rmpopt.merge_patch(input.at("rmp"));

    if (input.contains("molecule")) {
        // get the path of the system file
        std::filesystem::path syspath = inputpath / std::filesystem::path(input.at("molecule").at("file"));

        // assign the interface input path and define the result container
        orcaopt["folder"] = inputpath; Result res; Integrals ints(true);

        // make all the input and output paths absolute
        orcaopt.at("dynamics").at("output") = inputpath / orcaopt.at("dynamics").at("output");
        rhfopt.at("dynamics").at("output") = inputpath / rhfopt.at("dynamics").at("output");
        rmpopt.at("dynamics").at("output") = inputpath / rmpopt.at("dynamics").at("output");
        orcaopt.at("interface") = inputpath / orcaopt.at("interface");

        // create the system from the system file
        std::ifstream mstream(syspath); System system(mstream, input.at("molecule").at("basis")); mstream.close();

        // calculate all the atomic integrals
        if (!input.contains("orca")) {std::cout << std::endl;
            std::printf("OVERLAP INTEGRALS: "); auto stimer = Timer::Now(); ints.S = Integral::Overlap(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(stimer)).c_str());
            std::printf("KINETIC INTEGRALS: "); auto ttimer = Timer::Now(); ints.T = Integral::Kinetic(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(ttimer)).c_str());
            std::printf("NUCLEAR INTEGRALS: "); auto vtimer = Timer::Now(); ints.V = Integral::Nuclear(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(vtimer)).c_str());
            std::printf("COULOMB INTEGRALS: "); auto jtimer = Timer::Now(); ints.J = Integral::Coulomb(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(jtimer)).c_str());
        }

        // print and export the atomic integrals
        if (input.contains("integral")) {
            if (input.at("integral").contains("print")) {
                if (intopt.at("print").at("overlap")) std::cout << "\nOVERLAP INTEGRALS:\n" << ints.S << std::endl;
                if (intopt.at("print").at("kinetic")) std::cout << "\nKINETIC INTEGRALS:\n" << ints.T << std::endl;
                if (intopt.at("print").at("nuclear")) std::cout << "\nNUCLEAR INTEGRALS:\n" << ints.V << std::endl;
                if (intopt.at("print").at("coulomb")) std::cout << "\nCOULOMB INTEGRALS:\n" << ints.J << std::endl;
            }
            if (input.at("integral").contains("export")) {
                if (intopt.at("export").at("overlap")) EigenWrite(inputpath / "S.mat", ints.S);
                if (intopt.at("export").at("kinetic")) EigenWrite(inputpath / "T.mat", ints.T);
                if (intopt.at("export").at("nuclear")) EigenWrite(inputpath / "V.mat", ints.V);
                if (intopt.at("export").at("coulomb")) EigenWrite(inputpath / "J.mat", ints.J);
            }
        }

        if (input.contains("orca")) {
            if (input.at("orca").contains("dynamics")) {
                Orca(orcaopt).dynamics(system, ints, res);
            } else {
                throw std::runtime_error("YOU HAVE TO SPECIFY THE DYNAMICS BLOCK IN THE ORCA BLOCK");
            }

        } else if (input.contains("rhf")) {
            RestrictedHartreeFock rhf(rhfopt); res = rhf.run(system, ints);
            std::printf("\nRHF ENERGY: %.14f\n", res.Etot);

            if (input.at("rhf").contains("dynamics")) {
                rhf.dynamics(system, ints, res);
            } else if (input.at("rhf").contains("hessian")) {
                res = rhf.hessian(system, ints, res); std::cout << "\nRHF HESSIAN:\n" << res.rhf.H << "\n" << "RHF HESSIAN NORM: " << res.rhf.H.norm() << std::endl;
                std::cout << "\nRHF FREQUENCIES:\n" << Method<RestrictedHartreeFock>::frequency(system, res.rhf.H) << std::endl;
            } else if (input.at("rhf").contains("gradient")) {
                res = rhf.gradient(system, ints, res); std::cout << "\nRHF GRADIENT:\n" << res.rhf.G << "\n" << "RHF GRADIENT NORM: " << res.rhf.G.norm() << std::endl;
            }

            if (input.contains("rmp")) {
                RestrictedMollerPlesset rmp(rhfopt, rmpopt); ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C);
                res = rmp.run(system, ints, res); std::printf("\nRMP2 ENERGY: %.14f\n", res.Etot);

                if (input.at("rmp").contains("dynamics")) {
                    rmp.dynamics(system, ints, res);
                } else if (input.at("rmp").contains("hessian")) {
                    res = rmp.hessian(system, ints, res); std::cout << "\nRMP HESSIAN:\n" << res.rmp.H << "\n" << "RMP HESSIAN NORM: " << res.rmp.H.norm() << std::endl;
                    std::cout << "\nRMP FREQUENCIES:\n" << Method<RestrictedMollerPlesset>::frequency(system, res.rmp.H) << std::endl;
                } else if (input.at("rmp").contains("gradient")) {
                    res = rmp.gradient(system, ints, res); std::cout << "\nRMP2 GRADIENT:\n" << res.rmp.G << "\n" << "RMP2 GRADIENT NORM: " << res.rmp.G.norm() << std::endl;
                }
            }
        }
    } else if (input.contains("model")) {
        ModelSystem model(input.at("model").at("potential"), input.at("model").at("limits"), input.at("model").at("ngrid")); ModelSolver::Result res;

        if (input.contains("solve")) {
            ModelSolver msv(msvopt); res = msv.run(model);
            std::cout << "\nENERGIES:\n" << res.energy << std::endl;
            ModelSystem::SaveWavefunction("wavefunction.dat", res.r, res.states, res.energy);
            Matrix<> U(res.r.size(), 2); U << res.r, res.u; EigenWrite("U.mat", U);
        }
    }

    // finalize the integral engine and print the elapsed time
    std::printf("\nELAPSED TIME: %s\n", Timer::Format(Timer::Elapsed(programtimer)).c_str());
}

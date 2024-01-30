// include the methods
#include "restrictedconfigurationinteraction.h"

// model system and solver
#include "modelsolver.h"

// interfaces
#include "orca.h"

// include specific
#include "printer.h"
#include "default.h"

// include parsing and json libraries
#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>

// option structures loaders for RHF
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Hessian, step);

// option structures loaders for RCI
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Hessian, step);

// option structures loaders for RMP
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Hessian, step);

// option structures loaders for ORCA
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Dynamics, iters, step, output);

// option structures loaders for model method
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic::Dynamics, iters, step, output);

// option loaders
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic, dynamics, real, step, iters, nstate, thresh, optimize, guess);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options, dynamics, gradient, hessian, maxiter, thresh);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options, dynamics, gradient, hessian);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, dynamics, gradient, hessian, order);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic, dynamics, step, iters, guess);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options, dynamics, interface, method, folder);

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
    } nthread = program.get<int>("-n"); std::cout << std::fixed << std::setprecision(14); Result res;

    // get the path of the input
    auto inputpath = std::filesystem::current_path() / std::filesystem::path(program.get("input")).parent_path();

    // print program name and compiler version
    std::printf("QUANTUM ACORN (GCC %d.%d.%d, ", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);

    // print library versions
    std::printf("EIGEN %d.%d.%d, ", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
    std::printf("LIBINT %d.%d.%d)\n\n", libint2::major(), libint2::minor(), libint2::micro());

    // open the input file and parse the input
    std::ifstream istream(program.get("input")); nlohmann::json input = nlohmann::json::parse(istream); istream.close();

    // parse the default options
    auto rmpopt = nlohmann::json::parse(rmpoptstr), orcopt = nlohmann::json::parse(orcoptstr);
    auto intopt = nlohmann::json::parse(intoptstr), rhfopt = nlohmann::json::parse(rhfoptstr);
    auto msaopt = nlohmann::json::parse(msaoptstr), mdlopt = nlohmann::json::parse(mdloptstr);
    auto msnopt = nlohmann::json::parse(msnoptstr), rciopt = nlohmann::json::parse(rcioptstr);

    // assign the interface input path
    orcopt["folder"] = inputpath;

    // patch the input json file and apply defaults
    if (input.contains("integral")) intopt.merge_patch(input.at("integral"));
    if (input.contains("model")) mdlopt.merge_patch(input.at("model"));
    if (input.contains("orca")) orcopt.merge_patch(input.at("orca"));
    if (input.contains("rhf")) rhfopt.merge_patch(input.at("rhf"));
    if (input.contains("rci")) rciopt.merge_patch(input.at("rci"));
    if (input.contains("rmp")) rmpopt.merge_patch(input.at("rmp"));

    // patch the inputs for the modelsolver
    if (input.contains("solve") && mdlopt.at("potential").size() == 1) msaopt.merge_patch(input.at("solve"));
    if (input.contains("solve") && mdlopt.at("potential").size() != 1) msnopt.merge_patch(input.at("solve"));

    if (input.contains("molecule")) {
        // get the path of the system file
        std::filesystem::path syspath = inputpath / std::filesystem::path(input.at("molecule").at("file"));

        // define the integral container
        Integrals ints(true);

        // make all the input and output paths absolute
        orcopt.at("dynamics").at("output") = inputpath / orcopt.at("dynamics").at("output");
        rhfopt.at("dynamics").at("output") = inputpath / rhfopt.at("dynamics").at("output");
        rciopt.at("dynamics").at("output") = inputpath / rciopt.at("dynamics").at("output");
        rmpopt.at("dynamics").at("output") = inputpath / rmpopt.at("dynamics").at("output");
        orcopt.at("interface") = inputpath / orcopt.at("interface");

        // create the system from the system file
        std::ifstream mstream(syspath); System system(mstream, input.at("molecule").at("basis")); mstream.close();

        // print the integral title
        Printer::Title("INTEGRALS OVER ATOMIC ORBITALS");

        // calculate all the atomic integrals
        if (!input.contains("orca")) {
            std::printf("OVERLAP INTEGRALS: "); auto stimer = Timer::Now(); ints.S = Integral::Overlap(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(stimer)).c_str());
            std::printf("KINETIC INTEGRALS: "); auto ttimer = Timer::Now(); ints.T = Integral::Kinetic(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(ttimer)).c_str());
            std::printf("NUCLEAR INTEGRALS: "); auto vtimer = Timer::Now(); ints.V = Integral::Nuclear(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(vtimer)).c_str());
            std::printf("COULOMB INTEGRALS: "); auto jtimer = Timer::Now(); ints.J = Integral::Coulomb(system), std::printf("%s\n", Timer::Format(Timer::Elapsed(jtimer)).c_str());
        } std::cout << std::endl;

        // print and export the atomic integrals
        if (input.contains("integral")) {
            if (input.at("integral").contains("print")) {
                if (intopt.at("print").at("overlap")) Printer::Print(ints.S, "OVERLAP INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("kinetic")) Printer::Print(ints.T, "KINETIC INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("nuclear")) Printer::Print(ints.V, "NUCLEAR INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("coulomb")) Printer::Print(ints.J, "COULOMB INTEGRALS"), std::cout << "\n";
            } if (input.at("integral").contains("export")) {
                if (intopt.at("export").at("overlap")) EigenWrite(inputpath / "S.mat", ints.S);
                if (intopt.at("export").at("kinetic")) EigenWrite(inputpath / "T.mat", ints.T);
                if (intopt.at("export").at("nuclear")) EigenWrite(inputpath / "V.mat", ints.V);
                if (intopt.at("export").at("coulomb")) EigenWrite(inputpath / "J.mat", ints.J);
            }
        }

        // if the ORCA calculation is requested
        if (input.contains("orca")) {
            // if dynamics block is specified, run it else throw an error
            if (input.at("orca").contains("dynamics")) Orca(orcopt).dynamics(system, ints, res);
            else throw std::runtime_error("YOU HAVE TO DO DYNAMICS WITH ORCA");

        // if RHF calculation is requested
        } else if (input.contains("rhf")) {
            // print the title
            Printer::Title("RESTRICTED HARTREE-FOCK");

            // create the RHF object and run the calculation
            RestrictedHartreeFock rhf(rhfopt); res = rhf.run(system, ints);

            // print and export the RHF results
            if (input.at("rhf").contains("print")) {
                if (rhfopt.at("print").at("hcore")) Printer::Print(ints.T + ints.V, "CORE HAMILTONIAN MATRIX"), std::cout << "\n";
                if (rhfopt.at("print").at("coef")) Printer::Print(res.rhf.C, "COEFFICIENT MATRIX"), std::cout << "\n";
                if (rhfopt.at("print").at("density")) Printer::Print(res.rhf.D, "DENSITY MATRIX"), std::cout << "\n";
                if (rhfopt.at("print").at("orben")) Printer::Print(res.rhf.eps, "ORBITAL ENERGIES"), std::cout << "\n";
            } if (input.at("rhf").contains("export")) {
                if (rhfopt.at("export").at("hcore")) EigenWrite(inputpath / "H.mat", Matrix<>(ints.T + ints.V));
                if (rhfopt.at("export").at("coef")) EigenWrite(inputpath / "C.mat", res.rhf.C);
                if (rhfopt.at("export").at("density")) EigenWrite(inputpath / "D.mat", res.rhf.D);
                if (rhfopt.at("export").at("orben")) EigenWrite(inputpath / "EPS.mat", res.rhf.eps);
            }

            // print the total energy
            Printer::Print(res.Etot, "RESTRICTED HARTREE-FOCK ENERGY"), std::cout << "\n";

            // if the RHF dynamic block is used
            if (input.at("rhf").contains("dynamics")) {
                // run the dynamics from the created object and print newline
                rhf.dynamics(system, ints, res); std::cout << std::endl;

            // if the RHF hessian is requested
            } else if (input.at("rhf").contains("hessian")) {
                // calculate the hessian and print it
                res = rhf.hessian(system, ints, res); Printer::Print(res.rhf.H, "RESTRICTED HARTREE-FOCK HESSIAN"); std::cout << std::endl;

                // calculate the vibrational frequencies from the hessian and print them
                Printer::Print(Method<RestrictedHartreeFock>::frequency(system, res.rhf.H), "RESTRICTED HARTREE-FOCK FREQUENCIES"); std::cout << std::endl;

            // if the RHF gradient is requested
            } else if (input.at("rhf").contains("gradient")) {
                // calculate and print the gradient
                res = rhf.gradient(system, ints, res); Printer::Print(res.rhf.G, "RESTRICTED HARTREE-FOCK GRADIENT"); std::cout << std::endl;
            }

            // if configuratio interaction was requested
            if (input.contains("rci")) {
                // print the title
                Printer::Title("RESTRICTED CONFIGURATION INTERACTION");

                // create the RCI object and transform the integrals in the MS basis
                RestrictedConfigurationInteraction rci(rhfopt, rciopt); ints.Jms = Transform::CoulombSpin(ints.J, res.rhf.C);
                ints.Tms = Transform::SingleSpin(ints.T, res.rhf.C), ints.Vms = Transform::SingleSpin(ints.V, res.rhf.C);

                // perform the calculation
                res = rci.run(system, ints, res);

                // print and export the RCI prerequisities results
                if (input.at("rci").contains("print")) {
                    if (rciopt.at("print").at("kineticms")) Printer::Print(ints.Tms, "KINETIC INTEGRALS IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("nuclearms")) Printer::Print(ints.Vms, "NUCLEAR INTEGRALS IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("hcorems")) Printer::Print(ints.Tms + ints.Vms, "CORE HAMILTONIAN MATRIX IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("coulombms")) Printer::Print(ints.Jms, "COULOMB INTEGRALS IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("hamiltonian")) Printer::Print(res.rci.F, "CI HAMILTONIAN"), std::cout << "\n";
                    if (rciopt.at("print").at("energies")) Printer::Print(res.rci.eps, "EXCITED STATE ENERGIES"), std::cout << "\n";
                } if (input.at("rci").contains("export")) {
                    if (rciopt.at("export").at("kineticms")) EigenWrite(inputpath / "TMS.mat", ints.Tms);
                    if (rciopt.at("export").at("nuclearms")) EigenWrite(inputpath / "VMS.mat", ints.Vms);
                    if (rciopt.at("export").at("hcorems")) EigenWrite(inputpath / "HMS.mat", Matrix<>(ints.Tms + ints.Vms));
                    if (rciopt.at("export").at("coulombms")) EigenWrite(inputpath / "JMS.mat", ints.Jms);
                    if (rciopt.at("export").at("hamiltonian")) EigenWrite(inputpath / "HCI.mat", res.rci.F);
                    if (rciopt.at("export").at("energies")) EigenWrite(inputpath / "ECI.mat", res.rci.eps);
                }

                // print the resulting RCI energy
                Printer::Print(res.Etot, "RESTRICTED CONFIGURATION INTERACTION ENERGY"), std::cout << std::endl;

                // if the dynamics block is specified
                if (input.at("rci").contains("dynamics")) {
                    // perform the dynamics and print newline
                    rci.dynamics(system, ints, res); std::cout << std::endl;

                // if the hessian is requested
                } else if (input.at("rci").contains("hessian")) {
                    // calculate and print the hessian matrix
                    res = rci.hessian(system, ints, res); Printer::Print(res.rci.H, "RESTRICTED CONFIGURATION INTERACTION HESSIAN"); std::cout << std::endl;

                    // calculate and print the RCI vibrational frequencies
                    Printer::Print(Method<RestrictedConfigurationInteraction>::frequency(system, res.rci.H), "RESTRICTED CONFIGURATION INTERACTION FREQUENSIES"); std::cout << std::endl;

                // if gradient calculation was requested
                } else if (input.at("rci").contains("gradient")) {
                    // calculate and print the gradient matrix
                    res = rci.gradient(system, ints, res); Printer::Print(res.rci.G, "RESTRICTED CONFIGURATION INTERACTION GRADIENT"); std::cout << std::endl;
                }

            // if MP calculation was requested
            } else if (input.contains("rmp")) {
                // print the title
                Printer::Title("RESTRICTED MOLLER-PLESSET");

                // create the MP object and transform the coulomb tensor to MO basis
                RestrictedMollerPlesset rmp(rhfopt, rmpopt); ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C);

                // print and export the RMP prerequisities
                if (input.at("rmp").contains("print")) {
                    if (rmpopt.at("print").at("coulombmo")) Printer::Print(ints.Jmo, "COULOMB INTEGRALS IN MO BASIS"), std::cout << "\n";
                } if (input.at("rmp").contains("export")) {
                    if (rmpopt.at("export").at("coulombmo")) EigenWrite(inputpath / "JMO.mat", ints.Jmo);
                }

                // run the calculation and print the new line
                res = rmp.run(system, ints, res);

                // print the total energy
                Printer::Print(res.Etot, "RESTRICTED MOLLER-PLESSET ENERGY"); std::cout << std::endl;

                // if the MP dynamics was requested
                if (input.at("rmp").contains("dynamics")) {
                    // run the dynamics and print newline
                    rmp.dynamics(system, ints, res); std::cout << std::endl;

                // if the hessian calculation was requested
                } else if (input.at("rmp").contains("hessian")) {
                    // calculate the hessian and print the matrix
                    res = rmp.hessian(system, ints, res); Printer::Print(res.rmp.H, "RESTRICTED MOLLER-PLESSET HESSIAN"); std::cout << std::endl;

                    // calculate and print the RMP vibrational frequencies
                    Printer::Print(Method<RestrictedMollerPlesset>::frequency(system, res.rmp.H), "RESTRICTED MOLLER-PLESSET FREQUENSIES"); std::cout << std::endl;

                // if the RMP calculation was requested
                } else if (input.at("rmp").contains("gradient")) {
                    // perform the gradient calculation and print the matrix
                    res = rmp.gradient(system, ints, res); Printer::Print(res.rmp.G, "RESTRICTED MOLLER-PLESSET GRADIENT"); std::cout << std::endl;
                }
            }
        }

    // if the exact calculation was requested
    } else if (input.contains("model")) {
        // create the model system
        ModelSystem model(mdlopt.at("mass"), mdlopt.at("potential"), mdlopt.at("limits"), mdlopt.at("ngrid"));

        // if the solving was requested
        if (input.contains("solve")) {
            // print the title
            Printer::Title("EXACT QUANTUM DYNAMICS");

            // create the adianatic solver object
            ModelSolver msv(msaopt.get<ModelSolver::OptionsAdiabatic>());

            // if nonadiabatic model was specified assign to the solver the nonadiabatic version
            if (mdlopt.at("potential").size() > 1) msv = ModelSolver(msnopt.get<ModelSolver::OptionsNonadiabatic>());
            
            // run the calculation
            res = msv.run(model); 

            // if the optimization was performed print the resulting energies
            if (!msaopt.at("real") && mdlopt.at("potential").size() == 1) Printer::Print(res.msv.opten, "ENERGIES");

            // print new line
            std::cout << std::endl;

            // extract and save the potential points
            Matrix<> U(res.msv.r.size(), res.msv.U.cols() + 1); U << res.msv.r, res.msv.U; EigenWrite("U.mat", U);

        // throw error if nothing was specified for the model
        } else throw std::runtime_error("YOU HAVE TO DO SOMETHING WITH THE MODEL");
    }

    // finalize the integral engine and print the elapsed time
    std::stringstream ess; ess << "ELAPSED TIME: " << Timer::Format(Timer::Elapsed(programtimer)); Printer::Title(ess.str(), false);
}

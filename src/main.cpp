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

// code measuring define
#define MEASURE(T, F) std::cout << T << std::flush; {auto t = Timer::Now(); F; std::printf("%s\n", Timer::Format(Timer::Elapsed(t)).c_str());}

// option structures loaders for RHF
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Dynamics, iters, step, folder);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Hessian, step);

// option structures loaders for RCI
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Dynamics, iters, step, folder);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Hessian, step);

// option structures loaders for RMP
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Dynamics, iters, step, folder);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Hessian, step);

// option structures loaders for ORCA
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Dynamics, iters, step, folder);

// option structures loaders for model method
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic::Dynamics, iters, step, folder);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic::Dynamics, iters, step, folder);

// option loaders
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic, dynamics, real, step, iters, nstate, thresh, optimize, guess, folder);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options, dynamics, gradient, hessian, maxiter, thresh);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options, dynamics, gradient, hessian);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic, dynamics, step, iters, guess, folder);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, dynamics, gradient, hessian, order);
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
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
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

    // print the input file
    Printer::Title("INPUT FILE"); std::cout << input.dump(4) << std::endl << std::endl;

    // parse the default options
    auto intopt = nlohmann::json::parse(intoptstr);
    auto mdlopt = nlohmann::json::parse(mdloptstr);
    auto molopt = nlohmann::json::parse(moloptstr);
    auto msaopt = nlohmann::json::parse(msaoptstr);
    auto msnopt = nlohmann::json::parse(msnoptstr);
    auto orcopt = nlohmann::json::parse(orcoptstr);
    auto rciopt = nlohmann::json::parse(rcioptstr);
    auto rhfopt = nlohmann::json::parse(rhfoptstr);
    auto rmpopt = nlohmann::json::parse(rmpoptstr);

    // assign the dynamics output folder to the method options
    msaopt["dynamics"]["folder"] = inputpath;
    msnopt["dynamics"]["folder"] = inputpath;
    orcopt["dynamics"]["folder"] = inputpath;
    rciopt["dynamics"]["folder"] = inputpath;
    rhfopt["dynamics"]["folder"] = inputpath;
    rmpopt["dynamics"]["folder"] = inputpath;

    // assign the output paths to methods
    orcopt["folder"] = inputpath;
    msaopt["folder"] = inputpath;
    msnopt["folder"] = inputpath;

    // patch the input json file and apply defaults
    if (input.contains("integral")) intopt.merge_patch(input.at("integral"));
    if (input.contains("molecule")) molopt.merge_patch(input.at("molecule"));
    if (input.contains("model")) mdlopt.merge_patch(input.at("model"));
    if (input.contains("orca")) orcopt.merge_patch(input.at("orca"));
    if (input.contains("rci")) rciopt.merge_patch(input.at("rci"));
    if (input.contains("rhf")) rhfopt.merge_patch(input.at("rhf"));
    if (input.contains("rmp")) rmpopt.merge_patch(input.at("rmp"));

    // patch the inputs for the modelsolver
    if (input.contains("solve") && mdlopt.at("potential").size() == 1) msaopt.merge_patch(input.at("solve"));
    if (input.contains("solve") && mdlopt.at("potential").size() != 1) msnopt.merge_patch(input.at("solve"));

    // throw an error if restricted calculation cannot be performed due to multiplicity
    if (molopt.at("multiplicity") != 1 && input.contains("rci")) throw std::runtime_error("RESTRICTED CONFIGURATION INTERACTION CAN ONLY BE PERFORMED FOR SINGLET STATES");
    if (molopt.at("multiplicity") != 1 && input.contains("rmp")) throw std::runtime_error("RESTRICTED MOLLER-PLESSET CAN ONLY BE PERFORMED FOR SINGLET STATES");
    if (molopt.at("multiplicity") != 1 && input.contains("rhf")) throw std::runtime_error("RESTRICTED HARTREE-FOCK CAN ONLY BE PERFORMED FOR SINGLET STATES");

    // make all the interface paths absolute
    orcopt.at("interface") = inputpath / orcopt.at("interface");

    if (input.contains("molecule")) {
        // get the path of the system file
        std::filesystem::path syspath = inputpath / std::filesystem::path(input.at("molecule").at("file"));

        // define the integral container
        Integrals ints(true);

        // create the system from the system file
        std::ifstream mstream(syspath); System system(mstream, molopt.at("basis"), molopt.at("charge"), molopt.at("multiplicity")); mstream.close();

        // print the molecule specification
        Printer::Title("ATOM COORDINATE FILE AND BASIS"); std::cout << system << std::endl << std::endl;

        // calculate all the atomic integrals
        if (!input.contains("orca")) {Printer::Title("INTEGRALS OVER ATOMIC ORBITALS");
            MEASURE("OVERLAP INTEGRALS: ", ints.S = Integral::Overlap(system))
            MEASURE("KINETIC INTEGRALS: ", ints.T = Integral::Kinetic(system))
            MEASURE("NUCLEAR INTEGRALS: ", ints.V = Integral::Nuclear(system))
            MEASURE("COULOMB INTEGRALS: ", ints.J = Integral::Coulomb(system))
            if (input.contains("rhf") && input.at("rhf").contains("gradient")) {std::cout << std::endl;
                MEASURE("OVERLAP INTEGRAL DERIVATIVES: ", ints.dS = Integral::dOverlap(system))
                MEASURE("KINETIC INTEGRAL DERIVATIVES: ", ints.dT = Integral::dKinetic(system))
                MEASURE("NUCLEAR INTEGRAL DERIVATIVES: ", ints.dV = Integral::dNuclear(system))
                MEASURE("COULOMB INTEGRAL DERIVATIVES: ", ints.dJ = Integral::dCoulomb(system))
            } std::cout << std::endl;
        }

        // print and export the atomic integrals
        if (input.contains("integral")) {
            if (input.at("integral").contains("print")) {
                if (intopt.at("print").at("overlap")) Printer::Print(ints.S, "OVERLAP INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("kinetic")) Printer::Print(ints.T, "KINETIC INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("nuclear")) Printer::Print(ints.V, "NUCLEAR INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("coulomb")) Printer::Print(ints.J, "COULOMB INTEGRALS"), std::cout << "\n";
                if (intopt.at("print").at("doverlap")) Printer::Print(ints.dS, "OVERLAP INTEGRAL DERIVATIVES"), std::cout << "\n";
                if (intopt.at("print").at("dkinetic")) Printer::Print(ints.dT, "KINETIC INTEGRAL DERIVATIVES"), std::cout << "\n";
                if (intopt.at("print").at("dnuclear")) Printer::Print(ints.dV, "NUCLEAR INTEGRAL DERIVATIVES"), std::cout << "\n";
                if (intopt.at("print").at("dcoulomb")) Printer::Print(ints.dJ, "COULOMB INTEGRAL DERIVATIVES"), std::cout << "\n";
            } if (input.at("integral").contains("export")) {
                if (intopt.at("export").at("overlap")) EigenWrite(inputpath / "S.mat", ints.S);
                if (intopt.at("export").at("kinetic")) EigenWrite(inputpath / "T.mat", ints.T);
                if (intopt.at("export").at("nuclear")) EigenWrite(inputpath / "V.mat", ints.V);
                if (intopt.at("export").at("coulomb")) EigenWrite(inputpath / "J.mat", ints.J);
                if (intopt.at("export").at("doverlap")) EigenWrite(inputpath / "dS.mat", ints.dS);
                if (intopt.at("export").at("dkinetic")) EigenWrite(inputpath / "dT.mat", ints.dT);
                if (intopt.at("export").at("dnuclear")) EigenWrite(inputpath / "dV.mat", ints.dV);
                if (intopt.at("export").at("dcoulomb")) EigenWrite(inputpath / "dJ.mat", ints.dJ);
            }
        }

        // choose what calculation to run
        if (input.contains("orca")) {Printer::Title("ORCA DYNAMICS");
            // if dynamics block is specified, run it else throw an error
            if (input.at("orca").contains("dynamics")) {Orca(orcopt).dynamics(system, ints, res); std::cout << std::endl;}
            else throw std::runtime_error("YOU HAVE TO DO DYNAMICS WITH ORCA");

        } else if (input.contains("rhf")) {Printer::Title("RESTRICTED HARTREE-FOCK");
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
            if (input.at("rhf").contains("dynamics")) {Printer::Title("RESTRICTED HARTREE-FOCK DYNAMICS");
                rhf.dynamics(system, ints, res); std::cout << std::endl;

            // if the RHF hessian is requested
            } else if (input.at("rhf").contains("hessian")) {Printer::Title("RESTRICTED HARTREE-FOCK FREQUENCY CALCULATION");
                res = rhf.hessian(system, ints, res); Printer::Print(res.rhf.H, "RESTRICTED HARTREE-FOCK HESSIAN"); std::cout << std::endl;
                Printer::Print(rhf.frequency(system, res.rhf.H), "RESTRICTED HARTREE-FOCK FREQUENCIES"); std::cout << std::endl;

            // if the RHF gradient is requested
            } else if (input.at("rhf").contains("gradient")) {Printer::Title("RESTRICTED HARTREE-FOCK GRADIENT CALCULATION");
                res = rhf.gradient(system, ints, res); Printer::Print(res.rhf.G, "RESTRICTED HARTREE-FOCK GRADIENT"); std::cout << std::endl;
            }

            // if configuration interaction was requested
            if (input.contains("rci")) {Printer::Title("RESTRICTED CONFIGURATION INTERACTION");
                // create the RCI object
                RestrictedConfigurationInteraction rci(rhfopt, rciopt);

                // transform one-electron integrals to MS basis
                MEASURE("KINETIC INTEGRALS IN MS BASIS: ", ints.Tms = Transform::SingleSpin(ints.T, res.rhf.C))
                MEASURE("NUCLEAR INTEGRALS IN MS BASIS: ", ints.Vms = Transform::SingleSpin(ints.V, res.rhf.C))

                // transform the Coulomb integrals to MS basis
                MEASURE("COULOMB INTEGRALS IN MS BASIS: ", ints.Jms = Transform::CoulombSpin(ints.J, res.rhf.C))

                // newline
                std::cout << std::endl;

                // print and export the integrals in MS basis
                if (input.at("rci").contains("print")) {
                    if (rciopt.at("print").at("kineticms")) Printer::Print(ints.Tms, "KINETIC INTEGRALS IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("nuclearms")) Printer::Print(ints.Vms, "NUCLEAR INTEGRALS IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("hcorems")) Printer::Print(ints.Tms + ints.Vms, "CORE HAMILTONIAN MATRIX IN MS BASIS"), std::cout << "\n";
                    if (rciopt.at("print").at("coulombms")) Printer::Print(ints.Jms, "COULOMB INTEGRALS IN MS BASIS"), std::cout << "\n";
                } if (input.at("rci").contains("export")) {
                    if (rciopt.at("export").at("kineticms")) EigenWrite(inputpath / "TMS.mat", ints.Tms);
                    if (rciopt.at("export").at("nuclearms")) EigenWrite(inputpath / "VMS.mat", ints.Vms);
                    if (rciopt.at("export").at("hcorems")) EigenWrite(inputpath / "HMS.mat", Matrix<>(ints.Tms + ints.Vms));
                    if (rciopt.at("export").at("coulombms")) EigenWrite(inputpath / "JMS.mat", ints.Jms);
                }

                // perform the calculation
                res = rci.run(system, ints, res);

                // print and export the RCI results
                if (input.at("rci").contains("print")) {
                    if (rciopt.at("print").at("hamiltonian")) Printer::Print(res.rci.F, "CI HAMILTONIAN"), std::cout << "\n";
                    if (rciopt.at("print").at("energies")) Printer::Print(res.rci.eps, "EXCITED STATE ENERGIES"), std::cout << "\n";
                } if (input.at("rci").contains("export")) {
                    if (rciopt.at("export").at("hamiltonian")) EigenWrite(inputpath / "HCI.mat", res.rci.F);
                    if (rciopt.at("export").at("energies")) EigenWrite(inputpath / "ECI.mat", res.rci.eps);
                }

                // print the resulting RCI energy
                Printer::Print(res.Etot, "RESTRICTED CONFIGURATION INTERACTION ENERGY"), std::cout << std::endl;

                // if the dynamics block is specified
                if (input.at("rci").contains("dynamics")) {Printer::Title("RESTRICTED CONFIGURATION INTERACTION DYNAMICS");
                    rci.dynamics(system, ints, res); std::cout << std::endl;

                // if the hessian is requested
                } else if (input.at("rci").contains("hessian")) {Printer::Title("RESTRICTED CONFIGURATION FREQUENCY CALCULATION");
                    res = rci.hessian(system, ints, res); Printer::Print(res.rci.H, "RESTRICTED CONFIGURATION INTERACTION HESSIAN"); std::cout << std::endl;
                    Printer::Print(rci.frequency(system, res.rci.H), "RESTRICTED CONFIGURATION INTERACTION FREQUENCIES"); std::cout << std::endl;

                // if gradient calculation was requested
                } else if (input.at("rci").contains("gradient")) {Printer::Title("RESTRICTED CONFIGURATION INTERACTION GRADIENT CALCULATION");
                    res = rci.gradient(system, ints, res); Printer::Print(res.rci.G, "RESTRICTED CONFIGURATION INTERACTION GRADIENT"); std::cout << std::endl;
                }

            // if RMP calculation was requested
            } else if (input.contains("rmp")) {Printer::Title("RESTRICTED MOLLER-PLESSET");
                // create the RMP object
                RestrictedMollerPlesset rmp(rhfopt, rmpopt);

                // transform the coulomb tensor to MO basis
                MEASURE("COULOMB INTEGRALS IN MO BASIS: ", ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C)) std::cout << std::endl;

                // print and export the Coulomb integrals in MO
                if (input.at("rmp").contains("print")) {
                    if (rmpopt.at("print").at("coulombmo")) Printer::Print(ints.Jmo, "COULOMB INTEGRALS IN MO BASIS"), std::cout << "\n";
                } if (input.at("rmp").contains("export")) {
                    if (rmpopt.at("export").at("coulombmo")) EigenWrite(inputpath / "JMO.mat", ints.Jmo);
                }

                // run the calculation and print the new line
                res = rmp.run(system, ints, res);

                // print the total energy
                Printer::Print(res.Etot, "RESTRICTED MOLLER-PLESSET ENERGY"); std::cout << std::endl;

                // if the RMP dynamics was requested
                if (input.at("rmp").contains("dynamics")) {Printer::Title("RESTRICTED MOLLER-PLESSET DYNAMICS");
                    rmp.dynamics(system, ints, res); std::cout << std::endl;

                // if the hessian calculation was requested
                } else if (input.at("rmp").contains("hessian")) {Printer::Title("RESTRICTED MOLLER-PLESSET FREQUENCY CALCULATION");
                    res = rmp.hessian(system, ints, res); Printer::Print(res.rmp.H, "RESTRICTED MOLLER-PLESSET HESSIAN"); std::cout << std::endl;
                    Printer::Print(rmp.frequency(system, res.rmp.H), "RESTRICTED MOLLER-PLESSET FREQUENCIES"); std::cout << std::endl;

                // if the RMP calculation was requested
                } else if (input.at("rmp").contains("gradient")) {Printer::Title("RESTRICTED MOLLER-PLESSET GRADIENT CALCULATION");
                    res = rmp.gradient(system, ints, res); Printer::Print(res.rmp.G, "RESTRICTED MOLLER-PLESSET GRADIENT"); std::cout << std::endl;
                }
            }
        }

    // if the exact calculation was requested
    } else if (input.contains("model")) {
        // create the model system
        ModelSystem model(mdlopt.at("mass"), mdlopt.at("potential"), mdlopt.at("limits"), mdlopt.at("ngrid"));

        // if the solving was requested
        if (input.contains("solve")) {Printer::Title("EXACT QUANTUM DYNAMICS");

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
            Matrix<> U(res.msv.r.size(), res.msv.U.cols() + 1); U << res.msv.r, res.msv.U; EigenWrite(inputpath / "U.mat", U);
        }
    }

    // finalize the integral engine and print the elapsed time
    std::stringstream ess; ess << "ELAPSED TIME: " << Timer::Format(Timer::Elapsed(programtimer)); Printer::Title(ess.str(), false);
}

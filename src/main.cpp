// include the methods
#include "restrictedconfigurationinteraction.h"
#include "unrestrictedhartreefock.h"

// model solver and populations
#include "modelsolver.h"
#include "population.h"

// interfaces
#include "bagel.h"
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
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Dynamics, iters, step, berendsen);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options::Hessian, step);

// option structures loaders for RCI
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Dynamics, iters, step, berendsen);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options::Hessian, step);

// option structures loaders for RMP
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Dynamics, iters, step, berendsen);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options::Hessian, step);

// option structures loaders for UHF
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(UnrestrictedHartreeFock::Options::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(UnrestrictedHartreeFock::Options::Dynamics, iters, step, berendsen);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(UnrestrictedHartreeFock::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(UnrestrictedHartreeFock::Options::Hessian, step);

// option structures loaders for BAGEL
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Bagel::Options::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Bagel::Options::Dynamics, iters, step, berendsen);

// option structures loaders for ORCA
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options::Dynamics, iters, step, berendsen);

// option structures loaders for model methods
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic::Optimize, step, iters);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic::Spectrum, potential, window, normalize, zpesub);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic::Dynamics, iters, step, berendsen);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic::Dynamics::Berendsen, tau, temp, timeout);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic::Dynamics, iters, step, berendsen);

// option loaders
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsAdiabatic, dynamics, real, step, iters, nstate, optimize, guess, savewfn, spectrum);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedConfigurationInteraction::Options, dynamics, gradient, hessian, excitations, nstate, state);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(UnrestrictedHartreeFock::Options, dynamics, gradient, hessian, maxiter, thresh);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedHartreeFock::Options, dynamics, gradient, hessian, maxiter, thresh);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelSolver::OptionsNonadiabatic, dynamics, step, iters, guess, savewfn);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictedMollerPlesset::Options, dynamics, gradient, hessian, order);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Bagel::Options, dynamics, interface, nstate, state);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Orca::Options, dynamics, interface);

int main(int argc, char** argv) {
    // create the argument parser and the program timer
    argparse::ArgumentParser program("acorn", "0.2.0", argparse::default_arguments::none); auto programtimer = Timer::Now();

    // add the command line arguments arguments
    program.add_argument("input").help("-- Input file to specify the calculation.");
    program.add_argument("-d", "--defaults").help("-- Available options and defaults.").default_value(false).implicit_value(true);
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
    ip = std::filesystem::current_path() / std::filesystem::path(program.get("input")).parent_path();

    // print program name and compiler version
    std::printf("QUANTUM ACORN (GCC %d.%d.%d, ", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);

    // print library versions
    std::printf("EIGEN %d.%d.%d, ", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
    std::printf("LIBINT %d.%d.%d)\n\n", libint2::major(), libint2::minor(), libint2::micro());

    // open the input file and parse the input
    std::ifstream istream(program.get("input")); nlohmann::json input = nlohmann::json::parse(istream); istream.close();

    // print the input file and create the full input
    Printer::Title("INPUT FILE"); std::cout << input.dump(4) << std::endl << std::endl; nlohmann::json fullinput;

    // parse the default options
    auto intopt = nlohmann::json::parse(intoptstr);
    auto mdlopt = nlohmann::json::parse(mdloptstr);
    auto molopt = nlohmann::json::parse(moloptstr);
    auto msaopt = nlohmann::json::parse(msaoptstr);
    auto msnopt = nlohmann::json::parse(msnoptstr);
    auto orcopt = nlohmann::json::parse(orcoptstr);
    auto bglopt = nlohmann::json::parse(bgloptstr);
    auto rciopt = nlohmann::json::parse(rcioptstr);
    auto rhfopt = nlohmann::json::parse(rhfoptstr);
    auto rmpopt = nlohmann::json::parse(rmpoptstr);
    auto uhfopt = nlohmann::json::parse(uhfoptstr);

    // patch the input json file and apply defaults
    if (input.contains("integral")) intopt.merge_patch(input.at("integral")), fullinput["integral"] = intopt;
    if (input.contains("molecule")) molopt.merge_patch(input.at("molecule")), fullinput["molecule"] = molopt;
    if (input.contains("bagel")) bglopt.merge_patch(input.at("bagel")), fullinput["bagel"] = bglopt;
    if (input.contains("model")) mdlopt.merge_patch(input.at("model")), fullinput["model"] = mdlopt;
    if (input.contains("orca")) orcopt.merge_patch(input.at("orca")), fullinput["orca"] = orcopt;
    if (input.contains("rci")) rciopt.merge_patch(input.at("rci")), fullinput["rci"] = rciopt;
    if (input.contains("rhf")) rhfopt.merge_patch(input.at("rhf")), fullinput["rhf"] = rhfopt;
    if (input.contains("uhf")) uhfopt.merge_patch(input.at("uhf")), fullinput["uhf"] = uhfopt;
    if (input.contains("rmp")) rmpopt.merge_patch(input.at("rmp")), fullinput["rmp"] = rmpopt;

    // patch the inputs for the modelsolver
    if (input.contains("solve") && mdlopt.at("potential").size() == 1) {
        msaopt.merge_patch(input.at("solve"));
        if (msaopt.at("spectrum").at("potential").get<std::string>().empty()) {
            msaopt.at("spectrum").at("potential") = mdlopt.at("potential").at(0).at(0);
        } fullinput["solve"] = msaopt;
    } else if (input.contains("solve") && mdlopt.at("potential").size() != 1) msnopt.merge_patch(input.at("solve")), fullinput["solve"] = msnopt;

    // throw an error if restricted calculation cannot be performed due to multiplicity
    if (molopt.at("multiplicity") != 1 && input.contains("rci")) throw std::runtime_error("RESTRICTED CONFIGURATION INTERACTION CAN ONLY BE PERFORMED FOR SINGLET STATES");
    if (molopt.at("multiplicity") != 1 && input.contains("rmp")) throw std::runtime_error("RESTRICTED MOLLER-PLESSET CAN ONLY BE PERFORMED FOR SINGLET STATES");
    if (molopt.at("multiplicity") != 1 && input.contains("rhf")) throw std::runtime_error("RESTRICTED HARTREE-FOCK CAN ONLY BE PERFORMED FOR SINGLET STATES");

    // print the defaults if requested
    if (program.get<bool>("-d")) Printer::Title("DEFAULTS"), std::cout << fullinput.dump(4) << std::endl << std::endl;

    if (input.contains("molecule")) {
        // get the path of the system file and initialize integrals
        std::filesystem::path syspath = ip / std::filesystem::path(input.at("molecule").at("file")); Integrals ints(true);

        // create the system from the system file
        std::ifstream mstream(syspath); System system(mstream, molopt.at("basis"), molopt.at("charge"), molopt.at("multiplicity"));

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
                if (intopt.at("export").at("overlap")) EigenWrite(std::filesystem::path(ip) / "S.mat", ints.S);
                if (intopt.at("export").at("kinetic")) EigenWrite(std::filesystem::path(ip) / "T.mat", ints.T);
                if (intopt.at("export").at("nuclear")) EigenWrite(std::filesystem::path(ip) / "V.mat", ints.V);
                if (intopt.at("export").at("coulomb")) EigenWrite(std::filesystem::path(ip) / "J.mat", ints.J);
                if (intopt.at("export").at("doverlap")) EigenWrite(std::filesystem::path(ip) / "dS.mat", ints.dS);
                if (intopt.at("export").at("dkinetic")) EigenWrite(std::filesystem::path(ip) / "dT.mat", ints.dT);
                if (intopt.at("export").at("dnuclear")) EigenWrite(std::filesystem::path(ip) / "dV.mat", ints.dV);
                if (intopt.at("export").at("dcoulomb")) EigenWrite(std::filesystem::path(ip) / "dJ.mat", ints.dJ);
            }
        }

        // perform scans if movies provided
        if (input.contains("rmp") && mstream.peek() != EOF) {Printer::Title("RESTRICTED MOLLER-PLESSET SCAN");
            RestrictedMollerPlesset rmp(rhfopt, rmpopt); Matrix<> scan = rmp.scan(system, mstream, res); EigenWrite(ip / std::filesystem::path("scan.dat"), scan); std::cout << std::endl;
        } else if (input.contains("rci") && mstream.peek() != EOF) {Printer::Title("RESTRICTED CONFIGURATION INTERACTION SCAN");
            RestrictedConfigurationInteraction rci(rhfopt, rciopt); Matrix<> scan = rci.scan(system, mstream, res); EigenWrite(ip / std::filesystem::path("scan.dat"), scan); std::cout << std::endl;
        } else if (input.contains("rhf") && mstream.peek() != EOF) {Printer::Title("RESTRICTED HARTREE-FOCK SCAN");
            RestrictedHartreeFock rhf(rhfopt); Matrix<> scan = rhf.scan(system, mstream, res); EigenWrite(ip / std::filesystem::path("scan.dat"), scan); std::cout << std::endl;
        } else if (input.contains("uhf") && mstream.peek() != EOF) {Printer::Title("UNRESTRICTED HARTREE-FOCK SCAN");
            UnrestrictedHartreeFock uhf(uhfopt); Matrix<> scan = uhf.scan(system, mstream, res); EigenWrite(ip / std::filesystem::path("scan.dat"), scan); std::cout << std::endl;
        } else if (input.contains("bagel") && mstream.peek() != EOF) {Printer::Title("BAGEL SCAN");
            Bagel bagel(bglopt); Matrix<> scan = bagel.scan(system, mstream, res); EigenWrite(ip / std::filesystem::path("scan.dat"), scan); std::cout << std::endl;
        } else if (input.contains("orca") && mstream.peek() != EOF) {Printer::Title("ORCA SCAN");
            Orca orca(orcopt); Matrix<> scan = orca.scan(system, mstream, res); EigenWrite(ip / std::filesystem::path("scan.dat"), scan); std::cout << std::endl;

        // choose what calculation to run when not scanning
        } else if (input.contains("bagel")) {Printer::Title("BAGEL CALCULATION");
            // run the calculation
            Bagel bagel(bglopt); res = bagel.run(system, ints, res); std::cout << std::endl;

            // print the excited state energies
            Printer::Print(res.Eexc, "BAGEL EXCITED ENERGIES"); std::cout << std::endl;

            // print the total energy
            Printer::Print(res.Etot, "BAGEL ENERGY"), std::cout << std::endl;

            // if dynamics block is specified, run it
            if (input.at("bagel").contains("dynamics")) {Printer::Title("BAGEL DYNAMICS");
                bagel.dynamics(system, ints, res); std::cout << std::endl;

            // if the BAGEL hessian is requested
            } else if (input.at("bagel").contains("hessian")) {
                throw std::runtime_error("BAGEL HESSIAN CALCULATION NOT IMPLEMENTED");

            // if the BAGEL gradient is requested
            } else if (input.at("bagel").contains("gradient")) {Printer::Title("BAGEL GRADIENT CALCULATION");
                res = bagel.gradient(system, ints, res); std::cout << std::endl; Printer::Print(res.G, "BAGEL GRADIENT"); std::cout << std::endl;
            }

        } else if (input.contains("orca")) {Printer::Title("ORCA CALCULATION");
            // run the calculation
            Orca orca(orcopt); res = orca.run(system, ints, res); std::cout << std::endl;

            // print the total energy
            Printer::Print(res.Etot, "ORCA ENERGY"), std::cout << std::endl;

            // if dynamics block is specified, run it
            if (input.at("orca").contains("dynamics")) {Printer::Title("ORCA DYNAMICS");
                orca.dynamics(system, ints, res); std::cout << std::endl;

            // if the ORCA hessian is requested
            } else if (input.at("orca").contains("hessian")) {
                throw std::runtime_error("ORCA HESSIAN CALCULATION NOT IMPLEMENTED");

            // if the ORCA gradient is requested
            } else if (input.at("orca").contains("gradient")) {Printer::Title("ORCA GRADIENT CALCULATION");
                res = orca.gradient(system, ints, res); std::cout << std::endl; Printer::Print(res.G, "ORCA GRADIENT"); std::cout << std::endl;
            }

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
                if (rhfopt.at("export").at("hcore")) EigenWrite(std::filesystem::path(ip) / "H.mat", Matrix<>(ints.T + ints.V));
                if (rhfopt.at("export").at("coef")) EigenWrite(std::filesystem::path(ip) / "C.mat", res.rhf.C);
                if (rhfopt.at("export").at("density")) EigenWrite(std::filesystem::path(ip) / "D.mat", res.rhf.D);
                if (rhfopt.at("export").at("orben")) EigenWrite(std::filesystem::path(ip) / "EPS.mat", res.rhf.eps);
            }

            // print the mulliken charges
            if (rhfopt.at("mulliken")) Printer::Print(Population::Mulliken(system, ints, res.rhf.D), "MULLIKEN CHARGES"), std::cout << std::endl;

            // print the total energy
            Printer::Print(res.Etot, "RESTRICTED HARTREE-FOCK ENERGY"), std::cout << std::endl;

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
                    if (rciopt.at("export").at("kineticms")) EigenWrite(std::filesystem::path(ip) / "TMS.mat", ints.Tms);
                    if (rciopt.at("export").at("nuclearms")) EigenWrite(std::filesystem::path(ip) / "VMS.mat", ints.Vms);
                    if (rciopt.at("export").at("hcorems")) EigenWrite(std::filesystem::path(ip) / "HMS.mat", Matrix<>(ints.Tms + ints.Vms));
                    if (rciopt.at("export").at("coulombms")) EigenWrite(std::filesystem::path(ip) / "JMS.mat", ints.Jms);
                }

                // perform the calculation
                res = rci.run(system, ints, res);

                // print and export the RCI results
                if (input.at("rci").contains("print")) {
                    if (rciopt.at("print").at("hamiltonian")) Printer::Print(res.rci.F, "CI HAMILTONIAN"), std::cout << "\n";
                    if (rciopt.at("print").at("energies")) Printer::Print(res.rci.Eexc, "EXCITED STATE ENERGIES"), std::cout << "\n";
                } if (input.at("rci").contains("export")) {
                    if (rciopt.at("export").at("hamiltonian")) EigenWrite(std::filesystem::path(ip) / "HCI.mat", res.rci.F);
                    if (rciopt.at("export").at("energies")) EigenWrite(std::filesystem::path(ip) / "ECI.mat", res.rci.Eexc);
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
                    if (rmpopt.at("export").at("coulombmo")) EigenWrite(std::filesystem::path(ip) / "JMO.mat", ints.Jmo);
                }

                // run the calculation
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

        } else if (input.contains("uhf")) {Printer::Title("UNRESTRICTED HARTREE-FOCK");
            // create the RHF object and run the calculation
            UnrestrictedHartreeFock uhf(uhfopt); res = uhf.run(system, ints);

            // print and export the RHF results
            if (input.at("uhf").contains("print")) {
                if (uhfopt.at("print").at("hcore")) Printer::Print(ints.T + ints.V, "CORE HAMILTONIAN MATRIX"), std::cout << "\n";
                if (uhfopt.at("print").at("coefa")) Printer::Print(res.uhf.Ca, "ALPHA COEFFICIENT MATRIX"), std::cout << "\n";
                if (uhfopt.at("print").at("coefb")) Printer::Print(res.uhf.Cb, "BETA COEFFICIENT MATRIX"), std::cout << "\n";
                if (uhfopt.at("print").at("densitya")) Printer::Print(res.uhf.Da, "ALPHA DENSITY MATRIX"), std::cout << "\n";
                if (uhfopt.at("print").at("densityb")) Printer::Print(res.uhf.Db, "BETA DENSITY MATRIX"), std::cout << "\n";
                if (uhfopt.at("print").at("orbena")) Printer::Print(res.uhf.epsa, "ALPHA ORBITAL ENERGIES"), std::cout << "\n";
                if (uhfopt.at("print").at("orbenb")) Printer::Print(res.uhf.epsb, "BETA ORBITAL ENERGIES"), std::cout << "\n";
            } if (input.at("uhf").contains("export")) {
                if (uhfopt.at("export").at("hcore")) EigenWrite(std::filesystem::path(ip) / "H.mat", Matrix<>(ints.T + ints.V));
                if (uhfopt.at("export").at("coefa")) EigenWrite(std::filesystem::path(ip) / "CA.mat", res.uhf.Ca);
                if (uhfopt.at("export").at("coefb")) EigenWrite(std::filesystem::path(ip) / "CB.mat", res.uhf.Cb);
                if (uhfopt.at("export").at("densitya")) EigenWrite(std::filesystem::path(ip) / "DA.mat", res.uhf.Da);
                if (uhfopt.at("export").at("densityb")) EigenWrite(std::filesystem::path(ip) / "DB.mat", res.uhf.Db);
                if (uhfopt.at("export").at("orbena")) EigenWrite(std::filesystem::path(ip) / "EPSA.mat", res.uhf.epsa);
                if (uhfopt.at("export").at("orbenb")) EigenWrite(std::filesystem::path(ip) / "EPSB.mat", res.uhf.epsb);
            }

            // print the total energy
            Printer::Print(res.Etot, "UNRESTRICTED HARTREE-FOCK ENERGY"), std::cout << "\n";

            // if the UHF dynamic block is used
            if (input.at("uhf").contains("dynamics")) {Printer::Title("UNRESTRICTED HARTREE-FOCK DYNAMICS");
                uhf.dynamics(system, ints, res); std::cout << std::endl;

            // if the UHF hessian is requested
            } else if (input.at("uhf").contains("hessian")) {Printer::Title("UNRESTRICTED HARTREE-FOCK FREQUENCY CALCULATION");
                res = uhf.hessian(system, ints, res); Printer::Print(res.uhf.H, "UNRESTRICTED HARTREE-FOCK HESSIAN"); std::cout << std::endl;
                Printer::Print(uhf.frequency(system, res.uhf.H), "UNRESTRICTED HARTREE-FOCK FREQUENCIES"); std::cout << std::endl;

            // if the UHF gradient is requested
            } else if (input.at("uhf").contains("gradient")) {Printer::Title("UNRESTRICTED HARTREE-FOCK GRADIENT CALCULATION");
                res = uhf.gradient(system, ints, res); Printer::Print(res.uhf.G, "UNRESTRICTED HARTREE-FOCK GRADIENT"); std::cout << std::endl;
            }
        }

    } else if (input.contains("model")) {
        // create the model system
        ModelSystem model(mdlopt.at("mass"), mdlopt.at("potential"), mdlopt.at("limits"), mdlopt.at("ngrid"));

        if (input.contains("solve")) {Printer::Title("EXACT QUANTUM DYNAMICS");
            // create the adiabatic solver object
            ModelSolver msv(msaopt.get<ModelSolver::OptionsAdiabatic>());

            // if nonadiabatic model was specified assign to the solver the nonadiabatic version
            if (mdlopt.at("potential").size() > 1) msv = ModelSolver(msnopt.get<ModelSolver::OptionsNonadiabatic>());

            // run the optimization if requested
            if (mdlopt.at("potential").size() == 1 && input.at("solve").contains("optimize")) {
                auto opt = msaopt.get<ModelSolver::OptionsAdiabatic>(); opt.real = false;
                opt.iters = opt.optimize.iters, opt.step = opt.optimize.step;
                res = ModelSolver(opt).run(model);
            }
            
            // run the calculation
            res = msv.run(model, res);

            // if the optimization was performed print the resulting energies
            if (!msaopt.at("real") && mdlopt.at("potential").size() == 1) Printer::Print(res.msv.opten, "ENERGIES"), std::cout << std::endl;

            // print the spectral moments
            if (res.msv.spectra.size() && msaopt.at("spectrum").contains("moments")) {
                for (size_t i = 0; i < res.msv.spectra.size(); i++) {
                    double M0 = Numpy::Moment(res.msv.f, res.msv.spectra.at(i), 00, 0), M1 = Numpy::Moment(res.msv.f, res.msv.spectra.at(i), 00, 1);
                    if (msaopt.at("spectrum").at("moments") >= 0) std::printf("MOMENT #00 OF SPECTRUM #%d: % 2.2e\n", (int)i + 1, M0);
                    if (msaopt.at("spectrum").at("moments") >= 1) std::printf("MOMENT #01 OF SPECTRUM #%d: % 2.2e\n", (int)i + 1, M1);
                    for (int j = 2; j <= msaopt.at("spectrum").at("moments"); j++) {
                        double MJ = Numpy::Moment(res.msv.f, res.msv.spectra.at(i), M1, j); std::printf("MOMENT #%02d OF SPECTRUM #%d: % 2.2e\n", j, (int)i + 1, MJ);
                    }
                } std::cout << std::endl;
            }

            // extract and save the potential points
            if (StringContains(mdlopt.at("potential").at(0).at(0), 'y')) {
                Matrix<> U(res.msv.r.size() * res.msv.r.size(), res.msv.U.cols() + 2);
                for (int i = 0; i < res.msv.r.size(); i++) {
                    for (int j = 0; j < res.msv.r.size(); j++) {
                        U(i * res.msv.r.size() + j, 0) = res.msv.r(i), U(i * res.msv.r.size() + j, 1) = res.msv.r(j), U(i * res.msv.r.size() + j, 2) = res.msv.U(i * res.msv.r.size() + j);
                    }
                } EigenWrite(std::filesystem::path(ip) / "U.mat", U);
            } else {Matrix<> U(res.msv.r.size(), res.msv.U.cols() + 1); U << res.msv.r, res.msv.U; EigenWrite(std::filesystem::path(ip) / "U.mat", U);}

            // save the spectral analysis results if available
            if (res.msv.spectra.size()) {
                for (size_t i = 0; i < res.msv.spectra.size(); i++) {
                    Matrix<> A(res.msv.t.size(), 3); A << res.msv.t, res.msv.acfs.at(i).real(), res.msv.acfs.at(i).imag(); EigenWrite(std::filesystem::path(ip) / ("acf" + std::to_string(i) + ".mat"), A);
                    Matrix<> F(res.msv.f.size(), 2); F << res.msv.f, res.msv.spectra.at(i).cwiseAbs(); EigenWrite(std::filesystem::path(ip) / ("spectrum" + std::to_string(i) + ".mat"), F);
                }
            }
        }
    }

    // print the elapsed time
    std::stringstream ess; ess << "ELAPSED TIME: " << Timer::Format(Timer::Elapsed(programtimer)); Printer::Title(ess.str(), false);
}

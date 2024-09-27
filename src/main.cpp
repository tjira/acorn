#include                 "integral.h"
#include              "hartreefock.h"
#include            "mollerplesset.h"
#include           "coupledcluster.h"
#include          "quantumdynamics.h"
#include        "classicaldynamics.h"
#include "configurationinteraction.h"
#include               <argparse.hpp>
#include              <sys/utsname.h>

std::tuple<nlohmann::json, Input> parse_input(const std::filesystem::path& path) {
    // open the input file
    std::ifstream input_file_stream(path); if (!input_file_stream.good()) throw std::runtime_error("COULD NOT OPEN THE `" + path.string() + "` FILE");

    // parse the input file
    nlohmann::json input_json = nlohmann::json::parse(input_file_stream);

    // create the patched input
    nlohmann::json patched_input = default_input; patched_input.merge_patch(input_json);

    // return the input scruct
    return {input_json, patched_input.get<Input>()};
}

std::unique_ptr<argparse::ArgumentParser> parse_arguments(int argc, char** argv) {
    // extract the executable path and start the program timer
    std::filesystem::path executable_path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();

    // set the environment variable for the basis set location
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) setenv("LIBINT_DATA_PATH", executable_path.c_str(), true);

    // define the program
    std::unique_ptr<argparse::ArgumentParser> program = std::make_unique<argparse::ArgumentParser>("Acorn Quantum Package", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program->add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program->add_argument("-i", "--input").help("Input file to specify the calculations.").nargs(argparse::nargs_pattern::at_least_one);
    program->add_argument("-n", "--nthread").help("Number of threads to use during calculation.").default_value(1).scan<'i', int>();

    // parse the command line arguments
    try {program->parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program->get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program->get<bool>("-h")) {std::cout << program->help().str(); exit(EXIT_SUCCESS);}

    // assign the number of threads
    nthread = program->get<int>("-n");

    // return the program
    return program;
}

int main(int argc, char** argv) {
    // define the program timer and parse the command line arguments
    Timepoint program_timer = Timer::Now(); std::unique_ptr<argparse::ArgumentParser> program = parse_arguments(argc, argv);

    // print the program header
    std::printf("ACORN QUANTUM PACKAGE\n\n");

    // print the OS information
    struct utsname utsname; uname(&utsname); std::printf("PROGRAM EXECUTED ON: %s %s %s\n\n", utsname.sysname, utsname.release, utsname.version);

    // print the library versions
    std::printf("GCC: %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
    std::printf(", JSON: %d.%d.%d", NLOHMANN_JSON_VERSION_MAJOR, NLOHMANN_JSON_VERSION_MINOR, NLOHMANN_JSON_VERSION_PATCH);
    std::printf(", LIBINT: %s", LIBINT_VERSION);
    std::printf(", EIGEN: %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
    std::printf(", TORCH: %s\n\n", TORCH_VERSION);

    // print the compilation and execution timestamps
    std::printf("PROGRAM COMPILED: %s\nPROGRAM EXECUTED: %s\n", __TIMESTAMP__, Timer::Local().c_str());

    // loop over all input files
    for (const std::string& program_input : program->get<std::vector<std::string>>("-i")) {

        // print the input file being processed
        std::printf("\nPROCESSING INPUT FILE: %s\n\n", program_input.c_str());

        // parse the input file, initialize the system
        auto [input_json, input] = parse_input(program_input); System system; if (input_json.contains("system")) system = System(input.system);

        // define all the variables to store the results
        double energy_hf = 0, energy_mp = 0, energy_cc = 0, energy_ci = 0; torch::Tensor H_AO, S_AO, J_AO, H_MS, J_MS_AP, F_MS, C_MS, E_CI, C_CI;

        // extract all the method booleans
        bool do_hartree_fock = input_json.contains("hartree_fock");
        bool do_integral = input_json.contains("integral") || do_hartree_fock;
        bool do_configuration_interaction = do_hartree_fock && input_json.at("hartree_fock").contains("configuration_interaction");
        bool do_coupled_cluster = do_hartree_fock && input_json.at("hartree_fock").contains("coupled_cluster");
        bool do_moller_plesset = do_hartree_fock && input_json.at("hartree_fock").contains("moller_plesset");
        bool do_transform = do_hartree_fock && (do_configuration_interaction || do_coupled_cluster || do_moller_plesset);
        bool do_quantum_dynamics = input_json.contains("quantum_dynamics");
        bool do_classical_dynamics = input_json.contains("classical_dynamics");

        // print an error if the system multiplicity is greater than two or less than one
        if (input_json.contains("system") && (system.get_multi() < 1 || system.get_multi() > 2)) throw std::runtime_error("PROVIDED MULTIPLICITY NOT SUPPORTED");

        // print the initial system information
        if (input_json.contains("system")) printf("%s BASIS, %lu ATOMS, %d ELECTRONS, %d BASIS FUNCTIONS\n\n", system.get_basis().c_str(), system.get_atoms().size(), system.electrons(), system.basis_functions());

        // integral calculation
        if (do_integral) {

            // create the timer and print the calculation header
            Timepoint integral_calculation_timer = Timer::Now(); std::printf("INTEGRAL CALCULATION: "); std::flush(std::cout);

            // calculate the one- and two-electron integrals
            std::tie(H_AO, S_AO, J_AO) = Integral(input.integral).calculate(system.get_atoms(), system.get_shells());

            // print the time taken to calculate the integrals
            std::printf("%s\n", Timer::Format(Timer::Elapsed(integral_calculation_timer)).c_str());

            // wrte the integrals to the disk
            if (input.integral.data_export.hamiltonian || input.integral.data_export.coulomb || input.integral.data_export.overlap) {

                // create the export timer and print the header
                Timepoint integral_export_timer = Timer::Now(); std::printf("WRITING INTS TO DISK: "); std::flush(std::cout);

                // export the integrals
                if (input.integral.data_export.hamiltonian) Export::TorchTensorDouble("H_AO.mat", H_AO);
                if (input.integral.data_export.overlap    ) Export::TorchTensorDouble("S_AO.mat", S_AO);
                if (input.integral.data_export.coulomb    ) Export::TorchTensorDouble("J_AO.mat", J_AO);

                // print the time taken to export the integrals
                std::printf("%s\n", Timer::Format(Timer::Elapsed(integral_export_timer)).c_str());
            }

            // print the new line
            std::cout << std::endl;
        }

        // Hartree-Fock method
        if (do_hartree_fock) {

            // perform the calculation
            std::tie(F_MS, C_MS, energy_hf) = HartreeFock(input.hartree_fock).run(system, H_AO, S_AO, J_AO, torch::zeros_like(S_AO));

            // print the final energy
            std::printf("\nFINAL HF ENERGY: %.11f\n", energy_hf + system.nuclear_repulsion());
        }

        // integral transform
        if (do_transform) {

            // start the timer and print the header
            Timepoint transform_timer = Timer::Now(); std::printf("INTEGRAL TRANSFORM: "); std::flush(std::cout);

            // transform the integrals
            H_MS = Transform::SingleSpin(H_AO, C_MS), J_MS_AP = Transform::DoubleSpinAntsymPhys(J_AO, C_MS);

            // print the time taken to transform the integrals
            std::printf("%s\n", Timer::Format(Timer::Elapsed(transform_timer)).c_str());
        }

        // MP method
        if (do_moller_plesset) {

            // calculate the MP energy
            energy_mp = MollerPlesset(input.hartree_fock.moller_plesset).run(system, F_MS, J_MS_AP);

            // print the final MP energy
            std::printf("\nFINAL %s ENERGY: %.14f\n", MollerPlesset(input.hartree_fock.moller_plesset).get_name().c_str(), energy_hf + energy_mp + system.nuclear_repulsion());
        }

        // configuration interaction method
        if (do_configuration_interaction) {

            // calculate the CI energies and coefficients
            std::tie(E_CI, C_CI) = ConfigurationInteraction(input.hartree_fock.configuration_interaction).run(system, H_MS, J_MS_AP); energy_ci = E_CI.index({0}).item<double>() - energy_hf;

            // print the final CI energy
            std::printf("\nFINAL %s ENERGY: %.14f\n", ConfigurationInteraction(input.hartree_fock.configuration_interaction).get_name().c_str(), energy_hf + energy_ci + system.nuclear_repulsion());
        }

        // coupled cluster method
        if (do_coupled_cluster) {

            // calculate the CC energy
            energy_cc = CoupledCluster(input.hartree_fock.coupled_cluster).run(system, F_MS, J_MS_AP);

            // print the final CC energy
            std::printf("\nFINAL %s ENERGY: %.14f\n", CoupledCluster(input.hartree_fock.coupled_cluster).get_name().c_str(), energy_hf + energy_cc + system.nuclear_repulsion());
        }

        // perform the quantum dynamics calculation
        if (do_quantum_dynamics) QuantumDynamics(input.quantum_dynamics).run(input.wavefunction);

        // perform the classical dynamics calculation
        if (do_classical_dynamics) ClassicalDynamics(input.classical_dynamics).run(input.wavefunction);
    }

    // print the program timer
    std::printf("\nTOTAL EXECUTION TIME: %s\n", Timer::Format(Timer::Elapsed(program_timer)).c_str());
}

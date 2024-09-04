#include                 "integral.h"
#include              "hartreefock.h"
#include            "mollerplesset.h"
#include           "coupledcluster.h"
#include          "quantumdynamics.h"
#include "configurationinteraction.h"
#include               <argparse.hpp>

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
    program->add_argument("-i", "--input").help("Input file to specify the calculations.").default_value("input.json");

    // parse the command line arguments
    try {program->parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program->get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program->get<bool>("-h")) {std::cout << program->help().str(); exit(EXIT_SUCCESS);}

    // return the program
    return program;
}

int main(int argc, char** argv) {
    // define the program timer and parse the command line arguments
    Timepoint program_timer = Timer::Now(); std::unique_ptr<argparse::ArgumentParser> program = parse_arguments(argc, argv);

    // parse the input file, initialize the system
    auto [input_json, input] = parse_input(program->get("-i")); System system; if (input_json.contains("system")) system = System(input.system);

    // define all the variables to store the results
    double energy_hf{}, energy_mp{}, energy_cc{}, energy_ci{}; torch::Tensor H_AO, S_AO, J_AO, C_MO, H_MS, F_MS, J_MS_AP, E_CI, C_CI;

    // extract all the method booleans
    bool do_hartree_fock = input_json.contains("hartree_fock");
    bool do_integral = input_json.contains("integral") || do_hartree_fock;
    bool do_configuration_interaction = do_hartree_fock && input_json.at("hartree_fock").contains("configuration_interaction");
    bool do_coupled_cluster = do_hartree_fock && input_json.at("hartree_fock").contains("coupled_cluster");
    bool do_moller_plesset = do_hartree_fock && input_json.at("hartree_fock").contains("moller_plesset");
    bool do_transform = do_hartree_fock && (do_configuration_interaction || do_coupled_cluster || do_moller_plesset);
    bool do_quantum_dynamics = input_json.contains("quantum_dynamics");

    // print the initial system information
    if (input_json.contains("system")) printf("%s BASIS, %lu ATOMS, %d BASIS FUNCTIONS\n\n", system.get_basis().c_str(), system.get_atoms().size(), system.basis_functions());

    // calculate the integrals over the atomic orbitals
    if (do_integral) {
        std::tie(H_AO, S_AO, J_AO) = Integral(input.integral).calculate(system.get_atoms(), system.get_shells());
    }

    // perform the Hartree-Fock calculation to obtain the converged coefficients
    if (do_hartree_fock) {
        C_MO = HartreeFock(input.hartree_fock).run(system, H_AO, S_AO, J_AO, torch::zeros_like(S_AO));
        energy_hf = HartreeFock::get_energy(H_AO, HartreeFock::get_fock(H_AO, J_AO, HartreeFock::get_density(system, C_MO)), HartreeFock::get_density(system, C_MO));
        std::printf("\nFINAL HF ENERGY: %.14f\n", energy_hf + system.nuclear_repulsion());

    }

    // transform the Fock matrix, core Hamiltonian and Coulomb tensor to the basis of the molecular spinorbitals
    if (do_transform) {
        H_MS = Transform::SingleSpin(H_AO, C_MO);
        F_MS = Transform::SingleSpin(HartreeFock::get_fock(H_AO, J_AO, HartreeFock::get_density(system, C_MO)), C_MO);
        J_MS_AP = Transform::DoubleSpinAntsymPhys(J_AO, C_MO);
    }

    // calculate the MP2 correlation energy
    if (do_moller_plesset) {
        energy_mp = MollerPlesset(input.hartree_fock.moller_plesset).run(system, F_MS, J_MS_AP);
        std::printf("\nFINAL %s ENERGY: %.14f\n", MollerPlesset(input.hartree_fock.moller_plesset).get_name().c_str(), energy_hf + energy_mp + system.nuclear_repulsion());
    }

    // calculate the FCI energy
    if (do_configuration_interaction) {
        std::tie(E_CI, C_CI) = ConfigurationInteraction(input.hartree_fock.configuration_interaction).run(system, H_MS, J_MS_AP); energy_ci = E_CI.index({0}).item<double>() - energy_hf;
        std::printf("\nFINAL %s ENERGY: %.14f\n", ConfigurationInteraction(input.hartree_fock.configuration_interaction).get_name().c_str(), energy_hf + energy_ci + system.nuclear_repulsion());
    }

    // calculate the CCSD energy
    if (do_coupled_cluster) {
        energy_cc = CoupledCluster(input.hartree_fock.coupled_cluster).run(system, F_MS, J_MS_AP);
        std::printf("\nFINAL %s ENERGY: %.14f\n", CoupledCluster(input.hartree_fock.coupled_cluster).get_name().c_str(), energy_hf + energy_cc + system.nuclear_repulsion());
    }

    // perform the quantum dynamics calculation
    if (do_quantum_dynamics) {
        QuantumDynamics(input.quantum_dynamics).run(input.wavefunction);
    }

    // print the program timer
    std::printf("\nTOTAL EXECUTION TIME: %s\n", Timer::Format(Timer::Elapsed(program_timer)).c_str());
}

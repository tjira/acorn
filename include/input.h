#pragma once

#include <json.hpp>

struct Input {
    struct System {
        std::string basis, path;
    } system;
    struct Wavefunction {
        int dimension, grid_points; double mass, momentum;
        std::vector<double> grid_limits;
        std::vector<std::string> guess;
    } wavefunction;
    struct Integral {
        double precision;
    } integral;
    struct HartreeFock {
        int diis_size, max_iter;
        double threshold;
        struct ConfigurationInteraction {
            std::vector<int> excitation;
        } configuration_interaction;
        struct CoupledCluster {
            int max_iter;
            double threshold;
            std::vector<int> excitation, perturbation;
        } coupled_cluster;
        struct MollerPlesset {
            int order;
        } moller_plesset;
    } hartree_fock;
    struct QuantumDynamics {
        struct DataExport {
            bool diabatic_wavefunction, adiabatic_wavefunction, diabatic_density, adiabatic_density, diabatic_potential, adiabatic_potential, energy, position, momentum, acf;
        } data_export;
        bool align_wavefunction;
        int iterations, imaginary, real;
        double time_step;
        std::vector<std::vector<std::string>> potential;
    } quantum_dynamics;
};

extern nlohmann::json default_input;

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::System, basis, path);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::Wavefunction, dimension, mass, momentum, grid_limits, grid_points, guess);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::Integral, precision);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock::ConfigurationInteraction, excitation);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock::CoupledCluster, excitation, perturbation, max_iter, threshold);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock::MollerPlesset, order);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock, diis_size, max_iter, threshold, configuration_interaction, coupled_cluster, moller_plesset);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::QuantumDynamics::DataExport, diabatic_wavefunction, adiabatic_wavefunction, diabatic_density, adiabatic_density, diabatic_potential, adiabatic_potential, energy, position, momentum, acf);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::QuantumDynamics, potential, imaginary, real, iterations, time_step, data_export, align_wavefunction);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input, system, wavefunction, integral, hartree_fock, quantum_dynamics);

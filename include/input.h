#pragma once

#include <json.hpp>

struct Input {
    struct System {
        std::string basis, path; int charge, multiplicity;
    } system;
    struct Wavefunction {
        int dimension, grid_points; double mass, momentum;
        std::vector<double> grid_limits;
        std::vector<std::string> guess;
    } wavefunction;
    struct Integral {
        double precision;
        struct DataExport {
            bool hamiltonian, coulomb, overlap;
        } data_export;
    } integral;
    struct HartreeFock {
        bool generalized;
        int diis_size, max_iter;
        double threshold;
        struct ConfigurationInteraction {
            bool triplet;
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
            bool diabatic_wavefunction, adiabatic_wavefunction, diabatic_density, adiabatic_density, diabatic_population, adiabatic_population, diabatic_potential, adiabatic_potential, energy, position, momentum, acf;
        } data_export;
        int iterations, imaginary, real;
        double time_step;
        std::vector<std::vector<std::string>> potential;
    } quantum_dynamics;
    struct ClassicalDynamics {
        struct DataExport {
            bool diabatic_population, adiabatic_population, energy, energy_mean, position, position_mean, momentum, momentum_mean;
        } data_export;
        struct SurfaceHopping {
            int quantum_step_factor;
            std::string type;
        } surface_hopping;
        bool adiabatic;
        int iterations, trajectories, seed, log_interval_step, log_interval_traj;
        double time_step;
        std::vector<std::vector<std::string>> potential;
    } classical_dynamics;
};

extern nlohmann::json default_input; extern int nthread;

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::System, basis, path, charge, multiplicity);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::Wavefunction, dimension, mass, momentum, grid_limits, grid_points, guess);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::Integral::DataExport, hamiltonian, coulomb, overlap);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::Integral, precision, data_export);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock::ConfigurationInteraction, excitation, triplet);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock::CoupledCluster, excitation, perturbation, max_iter, threshold);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock::MollerPlesset, order);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::HartreeFock, diis_size, max_iter, threshold, configuration_interaction, coupled_cluster, moller_plesset, generalized);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::QuantumDynamics::DataExport, diabatic_wavefunction, adiabatic_wavefunction, diabatic_density, adiabatic_density, diabatic_population, adiabatic_population, diabatic_potential, adiabatic_potential, energy, position, momentum, acf);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::QuantumDynamics, potential, imaginary, real, iterations, time_step, data_export);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::ClassicalDynamics::DataExport, energy, position, momentum, diabatic_population, adiabatic_population, energy_mean, position_mean, momentum_mean);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::ClassicalDynamics::SurfaceHopping, type, quantum_step_factor);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input::ClassicalDynamics, potential, iterations, trajectories, time_step, data_export, seed, adiabatic, surface_hopping, log_interval_step, log_interval_traj);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Input, system, wavefunction, integral, hartree_fock, quantum_dynamics, classical_dynamics);

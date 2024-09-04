#pragma once

#include "export.h"
#include  "timer.h"

class QuantumDynamics {
public:
    struct IterationData {
        double energy, energy_error; std::complex<double> acf; Eigen::VectorXd position, momentum; Eigen::MatrixXd density_diabatic, density_adiabatic; Wavefunction diabatic_wavefunction, adiabatic_wavefunction;
    };

    QuantumDynamics(const Input::QuantumDynamics& input) : input(input) {}

    void export_trajectory(int iteration, const std::vector<IterationData>& iteration_data, const Eigen::MatrixXd& grid, bool imaginary) const;
    Eigen::MatrixXd get_diabatic_potential(const Eigen::MatrixXd& grid, const std::vector<std::string>& variables) const;
    std::tuple<std::vector<Eigen::MatrixXd>, Eigen::MatrixXd> get_transformation_matrices(const Eigen::MatrixXd& diabatic_potential) const;
    void print_iteration(int iteration, const IterationData& iteration_data, long elapsed) const;
    void run(const Input::Wavefunction& initial_diabatic_wavefunction_input) const;

private:
    Input::QuantumDynamics input;
};

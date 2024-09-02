#pragma once

#include "export.h"
#include  "timer.h"

class QuantumDynamics {
public:
    QuantumDynamics(const Input::QuantumDynamics& input) : input(input) {}

    std::vector<Eigen::MatrixXd> get_transformation_matrices(const Eigen::MatrixXd& diabatic_potential) const;
    Eigen::MatrixXd get_diabatic_potential(const Eigen::MatrixXd& grid, const std::vector<std::string>& variables) const;
    void run(const Wavefunction& initial_diabatic_wavefunction) const;

private:
    Input::QuantumDynamics input;
};

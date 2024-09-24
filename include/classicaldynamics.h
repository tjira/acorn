#pragma once

#include          "timer.h"
#include         "export.h"
#include    "landauzener.h"
#include   "wavefunction.h"
#include "fewestswitches.h"

class ClassicalDynamics {
public:
    struct TrajectoryData {
        Eigen::VectorXi state; Eigen::MatrixXd position, velocity; std::vector<Eigen::MatrixXd> diabatic_potential, adiabatic_potential;
    };

    ClassicalDynamics(const Input::ClassicalDynamics& input) : input(input) {}

    Eigen::MatrixXd calculate_derivative_coupling(const Eigen::MatrixXd& phi) const;
    Eigen::MatrixXd evaluate_potential(std::vector<std::vector<Expression>>& potential_expressions, const Eigen::VectorXd& position) const;
    Eigen::MatrixXd evaluate_potential_derivative(std::vector<std::vector<Expression>>& potential_expressions, const Eigen::VectorXd& position) const;
    void export_trajectories(const std::vector<TrajectoryData>& trajectory_data_vector, int mass) const;
    std::vector<std::vector<Expression>> get_potential_expression(const Wavefunction& initial_diabatic_wavefunction) const;
    void print_iteration(int trajectory, int iteration, const std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXi>& trajectory_data, double mass) const;
    void run(const Input::Wavefunction& initial_diabatic_wavefunction_input) const;

private:
    Input::ClassicalDynamics input;
};

namespace Export {
    void ClassicalTrajectories(const Input::ClassicalDynamics& input, const std::vector<ClassicalDynamics::TrajectoryData>& trajectory_data_vector, int mass);
}

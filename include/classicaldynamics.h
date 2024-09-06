#pragma once

#include       "export.h"
#include        "timer.h"
#include "wavefunction.h"
#include  "landauzener.h"
#include         <random>

class ClassicalDynamics {
public:
    struct TrajectoryData {
        Eigen::VectorXi state; Eigen::MatrixXd position, velocity; std::vector<Eigen::MatrixXd> diabatic_potential, adiabatic_potential;
    };

    ClassicalDynamics(const Input::ClassicalDynamics& input) : input(input) {}

    Eigen::MatrixXd evaluate_potential(std::vector<Expression>& potential_expressions, const Eigen::VectorXd& position) const;
    Eigen::MatrixXd evaluate_potential_derivative(std::vector<Expression>& potential_expressions, const Eigen::VectorXd& position) const;
    void export_trajectories(const std::vector<TrajectoryData>& trajectory_data_vector, int mass) const;
    std::vector<Expression> get_potential_expression(const Wavefunction& initial_diabatic_wavefunction) const;
    void print_iteration(int trajectory, int iteration, const std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXi>& trajectory_data, double mass) const;
    void run(const Input::Wavefunction& initial_diabatic_wavefunction_input) const;

private:
    Input::ClassicalDynamics input;
};

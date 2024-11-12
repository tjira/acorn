#pragma once

#include                           "input.h"
#include <unsupported/Eigen/MatrixFunctions>
#include                            <random>

// #include                         "export.h"

class FewestSwitches {
public:
    FewestSwitches(const Input::ClassicalDynamics::SurfaceHopping& input, bool adiabatic, int seed);

    static Eigen::MatrixXd calculate_derivative_coupling(const Eigen::MatrixXd& phi, const Eigen::MatrixXd& phi_prev, double time_step);
    static Eigen::MatrixXd calculate_derivative_coupling_kappa(const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, double time_step);
    static std::vector<std::tuple<int, double>> calculate_hopping_probabilities(const Eigen::VectorXcd& ci, const Eigen::MatrixXd& derivative_coupling, int state, double time_step);
    static Eigen::VectorXcd calculate_population_product(const Eigen::VectorXcd& population, const Eigen::VectorXd& potential, const Eigen::MatrixXd& derivative_coupling);
    static Eigen::VectorXcd propagate_population(const Eigen::VectorXcd& population, const Eigen::VectorXd& potential, const Eigen::MatrixXd& derivative_coupling, double time_step);

    std::tuple<Eigen::MatrixXd, Eigen::VectorXcd, int> jump(Eigen::VectorXcd population, const std::vector<Eigen::MatrixXd>& phi_vector, const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, int state, double time_step);

    // ~FewestSwitches() {
    //     Eigen::VectorXd data = Eigen::Map<Eigen::VectorXd>(offtdc.data(), offtdc.size());
    //
    //     Export::EigenMatrixDouble(std::string("OFFTDC-") + (input.kappa ? "K" : "") + "FSSH.mat", data, Eigen::VectorXd::LinSpaced(offtdc.size(), 1, offtdc.size()));
    // }
    // std::vector<double> offtdc;

private:
    bool adiabatic; std::uniform_real_distribution<double> dist; std::mt19937 mt; Input::ClassicalDynamics::SurfaceHopping input;
};

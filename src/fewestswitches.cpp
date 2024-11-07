#include "fewestswitches.h"

FewestSwitches::FewestSwitches(const Input::ClassicalDynamics::SurfaceHopping& input, bool adiabatic, int seed) : adiabatic(adiabatic), input(input) {
    // throw an error if diabatic mode requested
    if (!adiabatic && input.type == "fewest-switches") throw std::invalid_argument("DIABATIC MODE FOR FSSH NOT IMPLEMENTED YET");

    // initialize the random number generator
    this->dist = std::uniform_real_distribution<double>(0, 1), this->mt = std::mt19937(seed);
}

Eigen::MatrixXd FewestSwitches::calculate_derivative_coupling(const Eigen::MatrixXd& phi, const Eigen::MatrixXd& phi_prev, double time_step) {
    // define the derivative coupling matrix
    Eigen::MatrixXd derivative_coupling(phi.rows(), phi.cols());

    // calculate the overlap between the phi vectors
    for (int i = 0; i < phi.rows(); i++) {
        for (int j = 0; j < phi.cols(); j++) {
            derivative_coupling(i, j) = (phi_prev.col(i).transpose() * phi.col(j))(0);
        }
    }

    // return the derivative coupling matrix
    return derivative_coupling.log() / time_step;
}

Eigen::MatrixXd FewestSwitches::calculate_derivative_coupling_kappa(const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, double time_step) {
    // define the derivative coupling matrix
    Eigen::MatrixXd derivative_coupling(potential_vector.at(0).rows(), potential_vector.at(0).cols()); derivative_coupling.setZero();

    // calculate the derivative coupling
    for (int i = 0; i < derivative_coupling.rows(); i++) for (int j = 0; j < i; j++) {

        // calculate the energy differences needed for the second derivative
        double energy_difference_0 = potential_vector.at(iteration - 0)(i, i) - potential_vector.at(iteration - 0)(j, j);
        double energy_difference_1 = potential_vector.at(iteration - 1)(i, i) - potential_vector.at(iteration - 1)(j, j);
        double energy_difference_2 = potential_vector.at(iteration - 2)(i, i) - potential_vector.at(iteration - 2)(j, j);
        double energy_difference_3 = potential_vector.at(iteration - 3)(i, i) - potential_vector.at(iteration - 3)(j, j);

        // calculate the second derivative of the energy difference
        // double energy_difference_second_derivative = (energy_difference_0 - 2 * energy_difference_1 + energy_difference_2) / std::pow(time_step, 2);
        double energy_difference_second_derivative = (2 * energy_difference_0 - 5 * energy_difference_1 + 4 * energy_difference_2 - energy_difference_3) / std::pow(time_step, 2);

        // set the second derivative to zero if it is zero
        if (energy_difference_second_derivative < 0) energy_difference_second_derivative = 0;

        // calculate the derivative coupling
        derivative_coupling(i, j) = 0.5 * std::sqrt(energy_difference_second_derivative / energy_difference_0);
    }

    // return the derivative coupling matrix
    return derivative_coupling - derivative_coupling.transpose();
}

std::vector<std::tuple<int, double>> FewestSwitches::calculate_hopping_probabilities(const Eigen::VectorXcd& ci, const Eigen::MatrixXd& derivative_coupling, int state, double time_step) {
    // define the vector of hopping probabilities
    std::vector<std::tuple<int, double>> hopping_probabilities;

    // calculate the hopping probabilities
    for (int i = 0; i < ci.rows(); i++) {
        if (i != state) hopping_probabilities.push_back({i, -2 * (derivative_coupling(i, state) * ci(state) * ci.conjugate()(i)).real() / std::pow(std::abs(ci(state)), 2) * time_step});
    }

    // return the hopping probabilities
    return hopping_probabilities;
}

Eigen::VectorXcd FewestSwitches::calculate_population_product(const Eigen::VectorXcd& population, const Eigen::VectorXd& potential, const Eigen::MatrixXd& derivative_coupling) {
    return (potential.array() * population.array() / std::complex<double>(0, 1)).matrix() - derivative_coupling * population;
}

Eigen::VectorXcd FewestSwitches::propagate_population(const Eigen::VectorXcd& population, const Eigen::VectorXd& potential, const Eigen::MatrixXd& derivative_coupling, double time_step) {
    // calculate the intermediate Runge-Kutta steps
    Eigen::VectorXcd k1 = calculate_population_product(population,                      potential, derivative_coupling);
    Eigen::VectorXcd k2 = calculate_population_product(population + time_step * k1 / 2, potential, derivative_coupling);
    Eigen::VectorXcd k3 = calculate_population_product(population + time_step * k2 / 2, potential, derivative_coupling);
    Eigen::VectorXcd k4 = calculate_population_product(population + time_step * k3,     potential, derivative_coupling);

    // return the population
    return population + time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

std::tuple<Eigen::VectorXcd, int> FewestSwitches::jump(Eigen::VectorXcd population, const std::vector<Eigen::MatrixXd>& phi_vector, const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, int state, double time_step) {
    // define the derivative coupling and the new state
    Eigen::MatrixXd derivative_coupling; int new_state = state;

    // calculate the derivative coupling
    if (input.kappa && iteration > 2) {
        derivative_coupling = calculate_derivative_coupling_kappa(potential_vector, iteration, time_step);
    } else derivative_coupling = calculate_derivative_coupling(phi_vector.at(iteration), phi_vector.at(iteration - 1), time_step);

    // pripagate the populations
    for (int k = 0; k < input.quantum_step_factor; k++) {

        // propagate the population and generate a random number
        population = propagate_population(population, potential_vector.at(iteration).diagonal(), derivative_coupling, time_step / (double)input.quantum_step_factor); double random_number = dist(mt);

        // calculate the hopping probabilities
        std::vector<std::tuple<int, double>> hopping_probabilities = calculate_hopping_probabilities(population, derivative_coupling, state, time_step / (double)input.quantum_step_factor);

        // sort the transitions by probabilities from smallest to largest
        std::sort(hopping_probabilities.begin(), hopping_probabilities.end(), [](const auto& a, const auto& b) {return std::get<1>(a) < std::get<1>(b);});

        // add the cumulative probabilities
        for(size_t i = 0; i < hopping_probabilities.size(); i++) std::get<1>(hopping_probabilities.at(i)) += i ? std::get<1>(hopping_probabilities.at(i - 1)) : 0;

        // divide the probabilities by the maximum probability if the sum is greated than one
        for(size_t i = 0; i < hopping_probabilities.size() && std::get<1>(hopping_probabilities.back()) > 1; i++) std::get<1>(hopping_probabilities.at(i)) /= std::get<1>(hopping_probabilities.back());

        // check if a transition should be made
        for (auto& transition : hopping_probabilities) {
            if (new_state == state && random_number < std::get<1>(transition)) new_state = std::get<0>(transition);
        }
    }

    // return the propagated population and the new state
    return {population, new_state};
}

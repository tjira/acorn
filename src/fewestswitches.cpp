#include "fewestswitches.h"

FewestSwitches::FewestSwitches(const Input::ClassicalDynamics::SurfaceHopping& input, bool adiabatic, int seed) : input(input), adiabatic(adiabatic) {
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

std::tuple<Eigen::VectorXcd, int> FewestSwitches::jump(Eigen::VectorXcd population, const std::vector<Eigen::MatrixXd>& phi_vector, const Eigen::VectorXd& potential, int iteration, int state, double time_step) {
    // calculate the derivative coupling and define the new state
    Eigen::MatrixXd derivative_coupling = calculate_derivative_coupling(phi_vector.at(iteration), phi_vector.at(iteration - 1), time_step); int new_state = state;

    // pripagate the populations
    for (int k = 0; k < input.quantum_step_factor; k++) {

        // propagate the population and generate a random number
        population = propagate_population(population, potential, derivative_coupling, time_step / (double)input.quantum_step_factor); double random_number = dist(mt);

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

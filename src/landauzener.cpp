#include "landauzener.h"

int LandauZener::jump(const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, int state, double time_step, double random_number) const {
    // generate all the possible pairs of states and the transition container
    std::vector<std::vector<int>> pairs = Combinations(potential_vector.front().cols(), 2); std::vector<std::tuple<int, double, bool>> transitions;

    // loop over all the pairs of states
    for (size_t i = 0; i < pairs.size(); i++) {

        // define the current pair of states
        const int lower_state = pairs.at(i).at(0), upper_state = pairs.at(i).at(1);

        // skip if the current state is not in the combination and define the probability
        double transition_probability = 0; if (state != lower_state && state != upper_state) continue;

        // define the energy difference vector
        std::array<double, 3> energy_differences; std::array<double, 2> energy_differences_first_derivative; std::array<double, 1> energy_differences_second_derivative;

        // calculate the energy differences for the current and two previous iterations
        energy_differences.at(0) = potential_vector.at(iteration - 0)(upper_state, upper_state) - potential_vector.at(iteration - 0)(lower_state, lower_state);
        energy_differences.at(1) = potential_vector.at(iteration - 1)(upper_state, upper_state) - potential_vector.at(iteration - 1)(lower_state, lower_state);
        energy_differences.at(2) = potential_vector.at(iteration - 2)(upper_state, upper_state) - potential_vector.at(iteration - 2)(lower_state, lower_state);

        // calculate the first derivative of the energy difference
        energy_differences_first_derivative.at(0) = (energy_differences.at(0) - energy_differences.at(1)) / time_step;
        energy_differences_first_derivative.at(1) = (energy_differences.at(1) - energy_differences.at(2)) / time_step;

        // calculate the second derivative of the energy difference
        if (adiabatic) energy_differences_second_derivative.at(0) = (energy_differences.at(0) - 2 * energy_differences.at(1) + energy_differences.at(2)) / std::pow(time_step, 2);

        // calculate the adiabatic transition probaility
        if (adiabatic) transition_probability = std::exp(-0.5 * M_PI * std::sqrt(std::pow(energy_differences.at(0), 3) / energy_differences_second_derivative.at(0)));

        // calculate the diabatic transition probability
        if (!adiabatic) transition_probability = 1 - std::exp(-2 * M_PI * std::pow(potential_vector.at(iteration)(lower_state, upper_state), 2) / std::abs(energy_differences_first_derivative.at(0)));

        // add the transition to the container
        transitions.push_back(std::make_tuple(state == lower_state ? upper_state : lower_state, transition_probability, false));
 
        // check if the trajectory is at the jumping position for diabatic algorithm
        if (!adiabatic && energy_differences.at(0) * energy_differences.at(1) <= 0) std::get<2>(transitions.back()) = true;

        // check if the trajectory is at the jumping position for adiabatic algorithm
        if (adiabatic && energy_differences_first_derivative.at(0) * energy_differences_first_derivative.at(1) <= 0 && energy_differences_second_derivative.at(0) > 0) std::get<2>(transitions.back()) = true;
    }

    // sort the transitions by probabilities from smallest to largest
    std::sort(transitions.begin(), transitions.end(), [](const auto& a, const auto& b) {return std::get<1>(a) < std::get<1>(b);});

    // add the cumulative probabilities
    for(size_t i = 0; i < transitions.size(); i++) std::get<1>(transitions.at(i)) += i ? std::get<1>(transitions.at(i - 1)) : 0;

    // divide the probabilities by the maximum probability if the sum is greated than one
    for(size_t i = 0; i < transitions.size() && std::get<1>(transitions.back()) > 1; i++) std::get<1>(transitions.at(i)) /= std::get<1>(transitions.back());

    // return the new state if all conditions are met
    for (size_t i = 0; i < transitions.size(); i++) if (int new_state = std::get<0>(transitions.at(i)); std::get<2>(transitions.at(i))) {
        if (random_number > (i ? std::get<1>(transitions.at(i - 1)) : 0) && random_number < std::get<1>(transitions.at(i))) return new_state;
    }

    // return the same state
    return state;
}

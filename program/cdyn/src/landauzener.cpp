#include "landauzener.h"

LandauZener::LandauZener() = default;

inline std::vector<std::vector<int>> Combinations(int n, int k) {
    // create the bitmask that will get permuted and the resulting vector
    std::string bitmask(k, 1); bitmask.resize(n, 0); std::vector<std::vector<int>> combs;
 
    // generate the combinations
    do {std::vector<int> comb; comb.reserve(k);
        for (int j = 0; j < n; j++) {
            if (bitmask[j]) comb.push_back(j);
        } combs.push_back(comb);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    // return the result
    return combs;
}

LandauZener::LandauZener(int nstate, int points, bool adiabatic) : adiabatic(adiabatic) {
    combs = Combinations(nstate, 2); ed = Matrix::Zero(points, combs.size()), ded = Matrix::Zero(points, combs.size()), dded = Matrix::Zero(points, combs.size());;
}

std::vector<std::tuple<int, double, bool>> LandauZener::jump(const Matrix& U, int state, int i, double tstep) {
    // loop over all the state combinations
    for (size_t j = 0; j < combs.size(); j++) {

        // calculate the state energy differences
        ed(i, j) = U(combs.at(j).at(1), combs.at(j).at(1)) - U(combs.at(j).at(0), combs.at(j).at(0));

        // calculate the second derivative of the energy difference
        if (i > 1) dded(i, j) = (ed(i, j) - 2 * ed(i - 1, j) + ed(i - 2, j)) / (tstep * tstep);

        // calculate the first derivative of the energy difference
        if (i > 0) ded(i, j) = (ed(i, j) - ed(i - 1, j)) / tstep;
    }

    // define the container for the accepted transitions
    std::vector<std::tuple<int, double, bool>> transitions;

    // loop over all the state combinations
    for (size_t j = 0; j < combs.size() && i > 0; j++) {

        // skip if the current state is not in the combination and define the probability
        double P; if (state != combs.at(j).at(0) && state != combs.at(j).at(1)) continue;

        // calculate the probability of state change
        if (!adiabatic) P = 1 - std::exp(-2 * M_PI * std::pow(U(combs.at(j).at(0), combs.at(j).at(1)), 2) / std::abs(ded(i, j)));
        if ( adiabatic) P = std::exp(-0.5 * M_PI * std::sqrt(std::pow(ed(i, j), 3) / dded(i, j)));

        // add the transition to the container with the probability of the jump
        transitions.push_back(std::make_tuple(state == combs.at(j).at(0) ? combs.at(j).at(1) : combs.at(j).at(0), std::isnan(P) ? 0 : P, false));

        // if the trajectory is at a crossing point, enable the roll for the jump
        if (!adiabatic &&  ed(i, j) *  ed(i - 1, j) < 0) std::get<2>(transitions.back()) = true;
        if ( adiabatic && ded(i, j) * ded(i - 1, j) < 0) std::get<2>(transitions.back()) = true;
    }

    // sort the transitions by probabilities from smallest to largest
    std::sort(transitions.begin(), transitions.end(), [](const auto& a, const auto& b) {return std::get<1>(a) < std::get<1>(b);});

    // return transitions
    return transitions;
}

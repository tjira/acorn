#include "landauzener.h"
#include <iostream>

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
        // if (i > 2) dded(i, j) = (2 * ed(i, j) - 5 * ed(i - 1, j) + 4 * ed(i - 2, j) - ed(i - 3, j)) / (tstep * tstep);
        // if (i > 3) dded(i, j) = ((35.0/12) * ed(i, j) - (26.0/3) * ed(i - 1, j) + 9.5 * ed(i - 2, j) - (14.0/3) * ed(i - 3, j) + (11.0/12) * ed(i - 4, j)) / (tstep * tstep);

        // calculate the first derivative of the energy difference
        if (i > 0) ded(i, j) = (ed(i, j) - ed(i - 1, j)) / tstep;
        // if (i > 1) ded(i, j) = (1.5 * ed(i, j) - 2 * ed(i - 1, j) + 0.5 * ed(i - 2, j)) / tstep;
        // if (i > 2) ded(i, j) = ((11.0/6) * ed(i, j) - 3 * ed(i - 1, j) + 1.5 * ed(i - 2, j) - (1.0/3) * ed(i - 3, j)) / tstep;
    }

    // define the container for the accepted transitions
    std::vector<std::tuple<int, double, bool>> transitions;

    // loop over all the state combinations
    for (size_t j = 0; j < combs.size() && i > 0; j++) {

        // skip if the current state is not in the combination and define the probability
        double P = 0; if (state != combs.at(j).at(0) && state != combs.at(j).at(1)) continue;

        // calculate the probability of state change
        if (!adiabatic) {
            if (combs.size() == 3) {
                double l12 = std::sqrt(M_PI / std::abs(ded(i, 0))) * U(combs.at(0).at(0), combs.at(0).at(1));
                double l13 = std::sqrt(M_PI / std::abs(ded(i, 1))) * U(combs.at(1).at(0), combs.at(1).at(1));
                double l23 = std::sqrt(M_PI / std::abs(ded(i, 2))) * U(combs.at(2).at(0), combs.at(2).at(1));

                double P12 = 2 * std::pow(l12, 2) - 2 * l12 * l13 * l23 + std::pow(l13 * l23, 2) - std::pow(l12, 2) * (2 * std::pow(l12, 2) + 3 * std::pow(l13, 2) + 1 * std::pow(l23, 2));
                double P13 = 2 * std::pow(l13, 2) + 2 * l12 * l13 * l23 + std::pow(l12 * l23, 2) - std::pow(l13, 2) * (1 * std::pow(l12, 2) + 2 * std::pow(l13, 2) + 1 * std::pow(l23, 2));
                double P23 = 2 * std::pow(l23, 2) - 2 * l12 * l13 * l23 + std::pow(l12 * l13, 2) - std::pow(l23, 2) * (1 * std::pow(l12, 2) + 3 * std::pow(l13, 2) + 2 * std::pow(l23, 2));

                double P21 = 2 * std::pow(l12, 2) - 2 * l12 * l13 * l23 + std::pow(l13 * l23, 2) + std::pow(l12, 2) * (2 * std::pow(l12, 2) + 3 * std::pow(l13, 2) + 1 * std::pow(l23, 2));
                double P31 = 2 * std::pow(l13, 2) + 2 * l12 * l13 * l23 + std::pow(l12 * l23, 2) + std::pow(l13, 2) * (1 * std::pow(l12, 2) + 2 * std::pow(l13, 2) + 1 * std::pow(l23, 2));
                double P32 = 2 * std::pow(l23, 2) - 2 * l12 * l13 * l23 + std::pow(l12 * l13, 2) + std::pow(l23, 2) * (1 * std::pow(l12, 2) + 3 * std::pow(l13, 2) + 2 * std::pow(l23, 2));

                if (combs.at(j).at(0) == 0 && combs.at(j).at(1) == 1) P = state == combs.at(0).at(0) ? P12 : P21;
                if (combs.at(j).at(0) == 0 && combs.at(j).at(1) == 2) P = state == combs.at(1).at(0) ? P13 : P31;
                if (combs.at(j).at(0) == 1 && combs.at(j).at(1) == 2) P = state == combs.at(2).at(0) ? P23 : P32;
            } else {
                P = 1 - std::exp(-2 * M_PI * std::pow(U(combs.at(j).at(0), combs.at(j).at(1)), 2) / std::abs(ded(i, j)));
            }
        } else P = std::exp(-0.5 * M_PI * std::sqrt(std::pow(ed(i, j), 3) / dded(i, j)));

        // add the transition to the container with the probability of the jump
        transitions.push_back(std::make_tuple(state == combs.at(j).at(0) ? combs.at(j).at(1) : combs.at(j).at(0), std::isnan(P) ? 0 : P, false));

        // if the trajectory is at a crossing point, enable the roll for the jump
        if (!adiabatic &&  ed(i, j) *  ed(i - 1, j) <= 0                  ) std::get<2>(transitions.back()) = true;
        if ( adiabatic && ded(i, j) * ded(i - 1, j) <= 0 && dded(i, j) > 0) std::get<2>(transitions.back()) = true;
    }

    // sort the transitions by probabilities from smallest to largest
    std::sort(transitions.begin(), transitions.end(), [](const auto& a, const auto& b) {return std::get<1>(a) < std::get<1>(b);});

    // add the cumulative probabilities
    for(int i = 0; i < transitions.size(); i++) std::get<1>(transitions.at(i)) += i ? std::get<1>(transitions.at(i - 1)) : 0;

    // caluclate the sum of the probabilities
    double sum = 0; for(int i = 0; i < transitions.size(); i++) sum += std::get<1>(transitions.at(i));

    // divide the probabilities by the maximum probability if the sum is greated than one
    for(int i = 0; i < transitions.size() && sum > 1; i++) std::get<1>(transitions.at(i)) /= std::get<1>(transitions.back());

    // return transitions
    return transitions;
}

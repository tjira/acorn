#include "landauzener.h"

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

LandauZener::LandauZener(int nstate, bool adiabatic, int points, std::mt19937& mt) : adiabatic(adiabatic), mt(mt), dist(0, 1) {
    combs = Combinations(nstate, 2); de = Matrix::Zero(points, combs.size()), dz = Matrix::Zero(points, combs.size()), ddz = Matrix::Zero(points, combs.size());
}

int LandauZener::jump(const Matrix& U, int state, int i, double tstep) {
    // loop over all the state combinations
    for (size_t j = 0; j < combs.size(); j++) {

        // calculate the state energy differences
        de(i, j) = U(combs.at(j).at(1), combs.at(j).at(1)) - U(combs.at(j).at(0), combs.at(j).at(0));

        // calculate the second derivative
        if (i > 1) ddz(i, j) = (de(i, j) - 2 * de(i - 1, j) + de(i - 2, j)) / (tstep * tstep);
        // if (i > 2) ddz(i, j) = (2 * de(i, j) - 5 * de(i - 1, j) + 4 * de(i - 2, j) - de(i - 3, j)) / (tstep * tstep);
        // if (i > 3) ddz(i, j) = ((35.0/12) * de(i, j) - (26.0/3) * de(i - 1, j) + 9.5 * de(i - 2, j) - (14.0/3) * de(i - 3, j) + (11.0/12) * de(i - 4, j)) / (tstep * tstep);

        // calculate the first derivative
        if (i > 0) dz(i, j) = (de(i, j) - de(i - 1, j)) / tstep;
        // if (i > 1) dz(i, j) = (1.5 * de(i, j) - 2 * de(i - 1, j) + 0.5 * de(i - 2, j)) / tstep;
        // if (i > 2) dz(i, j) = ((11.0/6) * de(i, j) - 3 * de(i - 1, j) + 1.5 * de(i - 2, j) - (1.0/3) * de(i - 3, j)) / tstep;
    }

    // loop over all the state combinations
    for (size_t j = 0; j < combs.size() && i > 0; j++) {

        // skip if the current state is not in the combination and define the probability
        double P; if (state != combs.at(j).at(0) && state != combs.at(j).at(1)) continue;

        // calculate the probability of state change
        if (!adiabatic) P = 1 - std::exp(-2 * M_PI * std::pow(U(combs.at(j).at(0), combs.at(j).at(1)), 2) / std::abs(dz(i, j)));
        else P = std::exp(-0.5 * M_PI * std::sqrt(std::pow(de(i, j), 3) / ddz(i, j)));

        // return the new state if the jump is made
        if ((adiabatic && dz(i, j) * dz(i - 1, j) < 0 && dist(mt) < P) || (!adiabatic && de(i, j) * de(i - 1, j) < 0 && dist(mt) < P)) {
            return state == combs.at(j).at(0) ? combs.at(j).at(1) : combs.at(j).at(0);
        }
    }

    // return the current state if no jump is made
    return state;
}

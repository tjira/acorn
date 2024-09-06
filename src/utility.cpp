#include "utility.h"

std::vector<std::vector<int>> Combinations(int n, int k) {
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

std::vector<std::string> Split(const std::string& string, char delimeter) {
    std::vector<std::string> result; std::stringstream stringstream(string); std::string line; while (getline(stringstream, line, delimeter)) {result.push_back(line);} return result;
}

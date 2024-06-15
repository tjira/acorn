#pragma once

#include "linalg.h"

class LandauZener {
public:
    LandauZener(int nstate, bool adiabatic, int points, std::mt19937& mt); LandauZener() = default;

    // function to perform the Landau-Zener jump
    int jump(const Matrix& U, int state, int i, double tstep);

private:
    // general variables necessary for the Landau-Zener algorithm
    Matrix de, dz, ddz; std::vector<std::vector<int>> combs; bool adiabatic;

    // random number generator with uniform distribution
    std::mt19937 mt; std::uniform_real_distribution<double> dist;
};

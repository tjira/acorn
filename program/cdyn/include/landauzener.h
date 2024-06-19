#pragma once

#include "linalg.h"

class LandauZener {
public:
    LandauZener(int nstate, int points); LandauZener();

    // function to perform the Landau-Zener jump
    std::vector<std::tuple<int, double, bool>> jump(const Matrix& U, int state, int i, double tstep);

    // matrix getters
    const Matrix& getEd() const {return ed;} const Matrix& getDed() const {return ded;}

private:
    Matrix ed, ded; std::vector<std::vector<int>> combs;
};

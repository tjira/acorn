#pragma once

#include   "utility.h"
#include <Eigen/Dense>

class LandauZener {
public:
    LandauZener(bool adiabatic) : adiabatic(adiabatic) {}

    int jump(const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, int state, double time_step, double random_number) const;

private:
    bool adiabatic;
};

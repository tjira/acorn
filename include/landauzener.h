#pragma once

#include     "input.h"
#include   "utility.h"
#include <Eigen/Dense>
#include      <random>

class LandauZener {
public:
    LandauZener(const Input::ClassicalDynamics::SurfaceHopping& input, bool adiabatic, int seed);

    int jump(const std::vector<Eigen::MatrixXd>& potential_vector, int iteration, int state, double time_step);

private:
    bool adiabatic; std::uniform_real_distribution<double> dist; std::mt19937 mt; Input::ClassicalDynamics::SurfaceHopping input;
};

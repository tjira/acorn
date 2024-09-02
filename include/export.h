#pragma once

#include "wavefunction.h"

namespace Export {
    void EigenMatrixDouble(const std::string& path, const Eigen::MatrixXd& matrix, const Eigen::VectorXd& independent_variable = Eigen::VectorXd::Zero(0));
    void WavefunctionTrajectory(const std::string& path, const std::vector<Wavefunction>& wavefunction_trajectory, const Eigen::VectorXd& independent_variable = Eigen::VectorXd::Zero(0));
}

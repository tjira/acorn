#pragma once

#include "wavefunction.h"

namespace Export {
    void EigenMatrixDouble(const std::string& path, const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& independent_variable = Eigen::VectorXd::Zero(0));
}

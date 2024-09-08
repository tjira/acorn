#pragma once

#include <torch/torch.h>
#include   <Eigen/Dense>
#include       <fstream>
#include       <iomanip>

namespace Export {
    void EigenMatrixDouble(const std::string& path, const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& independent_variable = Eigen::VectorXd::Zero(0));
    void TorchTensorDouble(const std::string& path, const torch::Tensor& tensor);
}

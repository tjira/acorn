#pragma once

#include    <exprtk.hpp>
#include <torch/torch.h>
#include   <Eigen/Dense>

class Expression {
public:
    Expression() = default; Expression(const std::string& expression_string, const std::vector<std::string>& variable_strings);

    Eigen::VectorXd evaluate(const Eigen::MatrixXd& r);
    torch::Tensor evaluate(const torch::Tensor& r);

private:
    std::string expression_string; exprtk::expression<double> expression; std::vector<double> variables; std::vector<std::string> variable_strings;
};

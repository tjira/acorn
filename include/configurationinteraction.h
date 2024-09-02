#pragma once

#include "system.h"
#include  "timer.h"

class ConfigurationInteraction {
public:
    ConfigurationInteraction(const Input::HartreeFock::ConfigurationInteraction& input) : input(input) {}

    static std::tuple<std::vector<int>, std::vector<int>, int> align_determinants(std::vector<int> deta, const std::vector<int>& detb);
    static std::vector<std::vector<int>> generate_all_configurations(const System& system);
    static std::tuple<std::vector<int>, std::vector<int>> get_common_and_unique_spinorbitals(std::vector<int> deta, const std::vector<int>& detb);
    static double slater_condon_rules(std::vector<int> deta, const std::vector<int>& detb, const torch::TensorAccessor<double, 2>& H_MS_accessor, const torch::TensorAccessor<double, 4>& J_MS_AP_accessor);

    std::string get_name() const;
    std::tuple<torch::Tensor, torch::Tensor> run(const System& system, const torch::Tensor& H_MS, const torch::Tensor& J_MS_AP) const;

private:
    Input::HartreeFock::ConfigurationInteraction input;
};

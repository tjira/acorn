#pragma once

#include    "system.h"
#include     "timer.h"
#include "transform.h"
#include   "utility.h"

class MollerPlesset {
public:
    MollerPlesset(const Input::HartreeFock::MollerPlesset& input) : input(input) {}

    static double evaluate_contraction(const System& system, const std::string& contraction_string, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP, int order);

    std::string get_name() const;
    double run(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const;

private:
    Input::HartreeFock::MollerPlesset input;
};

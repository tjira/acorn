#pragma once

#include "system.h"
#include  "timer.h"

class HartreeFock {
public:
    HartreeFock(const Input::HartreeFock& input) : input(input) {}

    static torch::Tensor get_density(const System& system, const torch::Tensor& C_MO);
    static double get_energy(const torch::Tensor& H_AO, const torch::Tensor& F_AO, const torch::Tensor& D_MO);
    static torch::Tensor get_fock(const torch::Tensor& H_AO, const torch::Tensor& J_AO, const torch::Tensor& D_MO);

    torch::Tensor run(const System& system, const torch::Tensor& H_AO, const torch::Tensor& S_AO, const torch::Tensor& J_AO, torch::Tensor Dmo) const;

private:
    Input::HartreeFock input;
};

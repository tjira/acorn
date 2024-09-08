#pragma once

#include    "system.h"
#include     "timer.h"
#include "transform.h"

class HartreeFock {
public:
    HartreeFock(const Input::HartreeFock& input) : input(input) {}

    torch::Tensor get_density(const System& system, const torch::Tensor& C_MO) const;
    double get_energy(const torch::Tensor& H_AO, const torch::Tensor& F_AO, const torch::Tensor& D_MO) const;
    torch::Tensor get_fock(const torch::Tensor& H_AO, const torch::Tensor& J_AO, const torch::Tensor& D_MO) const;
    std::tuple<torch::Tensor, torch::Tensor, double> run(const System& system, const torch::Tensor& H_AO, const torch::Tensor& S_AO, const torch::Tensor& J_AO, torch::Tensor D_MO) const;
    std::tuple<torch::Tensor, torch::Tensor, double> scf(const System& system, const torch::Tensor& H, const torch::Tensor& S, const torch::Tensor& J, torch::Tensor D) const;

private:
    Input::HartreeFock input;
};

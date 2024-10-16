#pragma once

#include       "input.h"
#include   <libint2.hpp>
#include <torch/torch.h>

class Integral {
public:
    Integral(const Input::Integral& input) : precision(input.precision) {}

    // basic integrals
    static torch::Tensor double_electron(libint2::Engine& engine, const libint2::BasisSet& shells);
    static torch::Tensor single_electron(libint2::Engine& engine, const libint2::BasisSet& shells);

    // first derivatives
    static torch::Tensor double_electron_d1(libint2::Engine& engine, const libint2::BasisSet& shells, const std::vector<libint2::Atom>& atoms);
    static torch::Tensor single_electron_d1(libint2::Engine& engine, const libint2::BasisSet& shells, const std::vector<libint2::Atom>& atoms);

    // calculate the basic integrals
    std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor> calculate(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) const;

    // calculate the first derivatives
    std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor> calculate_d1(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) const;

private:
    double precision;
};

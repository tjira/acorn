#pragma once

#include       "input.h"
#include   <libint2.hpp>
#include <torch/torch.h>

class Integral {
public:
    Integral(const Input::Integral& input) : precision(input.precision) {}

    static torch::Tensor double_electron(libint2::Engine& engine, const libint2::BasisSet& shells);
    static torch::Tensor single_electron(libint2::Engine& engine, const libint2::BasisSet& shells);

    std::tuple<torch::Tensor, torch::Tensor, torch::Tensor> calculate(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) const;

private:
    double precision;
};

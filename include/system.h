#pragma once

#include       "input.h"
#include    "constant.h"
#include   <libint2.hpp>
#include <torch/torch.h>

using namespace torch::indexing;

class System {
public:
    System() = default; System(const Input::System& input);

    int basis_functions() const;
    int electrons() const;
    std::vector<libint2::Atom> get_atoms() const;
    std::string get_basis() const;
    int get_multi() const;
    libint2::BasisSet get_shells() const;
    double nuclear_repulsion() const;
    int virtual_spinorbitals() const;

private:
    int charge, multi; std::string basis; std::vector<int> atomic_numbers; torch::Tensor coordinates;
};

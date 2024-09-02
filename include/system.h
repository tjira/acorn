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
    std::vector<libint2::Atom> get_atoms() const;
    std::string get_basis() const;
    libint2::BasisSet get_shells() const;
    double nuclear_repulsion() const;
    int occupied_spatial_orbitals() const;
    int occupied_spinorbitals() const;
    int virtual_spinorbitals() const;

private:
    std::string basis; std::vector<int> atomic_numbers; torch::Tensor coordinates;
};

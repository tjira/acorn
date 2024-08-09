#pragma once

#include "tensor.h"
#include <libint2.hpp>

namespace Integral {
    // general integrals of the form (i|O|j) and (ij|O|kl) in chemists' notation
    torch::Tensor Single(libint2::Engine& engine, const libint2::BasisSet& shells);
    torch::Tensor Double(libint2::Engine& engine, const libint2::BasisSet& shells);

    // nuclear integral with the operator dependent of system coordinates
    torch::Tensor Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells);

    // kinetic, overlap, and coulomb integrals
    torch::Tensor Kinetic(const libint2::BasisSet& shells);
    torch::Tensor Overlap(const libint2::BasisSet& shells);
    torch::Tensor Coulomb(const libint2::BasisSet& shells);
}

#pragma once

#include "linalg.h"
#include <libint2.hpp>

namespace Integral {
    // general integrals of the form (i|O|j) and (ij|O|kl) in chemist's notation
    Matrix    Single(libint2::Engine& engine, const libint2::BasisSet& shells);
    Tensor<4> Double(libint2::Engine& engine, const libint2::BasisSet& shells);

    // nuclear integral with the operator dependent of system coordinates
    Matrix Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells);

    // kinetic, overlap, and coulomb integrals
    Matrix    Kinetic(const libint2::BasisSet& shells);
    Matrix    Overlap(const libint2::BasisSet& shells);
    Tensor<4> Coulomb(const libint2::BasisSet& shells);
}

#pragma once

#include "system.h"
#include "ptable.h"

#include <libint2.hpp>

namespace Integral {
    // general integrals of the form (i|O|j) and (ij|O|kl) in chemist's notation
    EigenMatrix<> Single(libint2::Engine& engine, const libint2::BasisSet& shells);
    EigenTensor<> Double(libint2::Engine& engine, const libint2::BasisSet& shells);

    // nuclear integral with the operator dependent of system coordinates
    EigenMatrix<> Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells);

    // kinetic, overlap, and coulomb integrals
    EigenMatrix<> Kinetic(const libint2::BasisSet& shells);
    EigenMatrix<> Overlap(const libint2::BasisSet& shells);
    EigenTensor<> Coulomb(const libint2::BasisSet& shells);
}

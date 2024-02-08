#pragma once

#include "system.h"

namespace libint2 {
    class Engine;
}

namespace Integral {
    // integrals over atomic orbitals
    Tensor<> Coulomb(const System& system);
    Matrix<> Kinetic(const System& system);
    Matrix<> Overlap(const System& system);
    Matrix<> Nuclear(const System& system);

    // integral derivatives over atomic orbitals
    Tensor<5> dCoulomb(const System& system);
    Tensor<3> dKinetic(const System& system);
    Tensor<3> dOverlap(const System& system);
    Tensor<3> dNuclear(const System& system);

    // general integrals
    Matrix<> Single(libint2::Engine& engine, const System& system);
    Tensor<> Double(libint2::Engine& engine, const System& system);

    // general integral derivatives
    Tensor<5> dDouble(libint2::Engine& engine, const System& system);
    Tensor<3> dSingle(libint2::Engine& engine, const System& system);

// integral container
} struct Integrals {Integrals(bool main = false), ~Integrals(); bool main; Matrix<> T, S, V, Tmo, Smo, Vmo, Tms, Sms, Vms; Tensor<3> dT, dS, dV; Tensor<4> J, Jmo, Jms; Tensor<5> dJ;};

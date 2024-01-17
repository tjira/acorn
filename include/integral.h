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

    // general integrals
    Matrix<> Single(libint2::Engine& engine, const System& system);
    Tensor<> Double(libint2::Engine& engine, const System& system);

    // initializer and finalizer
    void Initialize(); void Finalize();

// integral container
} struct Integrals {Matrix<> T, S, V; Tensor<> J, Jmo;};

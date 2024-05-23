#pragma once

#include "eigen.h"

namespace Transform {
    // two electron integrals
    EigenTensor<> CoulombSpin(const EigenTensor<>& Jao, const EigenMatrix<>& Cmo);
    EigenTensor<> CoulombSpatial(const EigenTensor<>& Jao, EigenMatrix<>& Cmo);

    // one electron integrals
    EigenMatrix<> SingleSpin(const EigenMatrix<>& Aao, const EigenMatrix<>& Cmo);
    EigenMatrix<> SingleSpatial(const EigenMatrix<>& Aao, EigenMatrix<>& Cmo);
}

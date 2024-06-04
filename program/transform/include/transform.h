#pragma once

#include "linalg.h"

namespace Transform {
    // two electron integrals
    Tensor<4> CoulombSpin(const Tensor<4>& Jao, const Matrix& Cmo);
    Tensor<4> CoulombSpatial(const Tensor<4>& Jao, Matrix& Cmo);

    // one electron integrals
    Matrix SingleSpin(const Matrix& Aao, const Matrix& Cmo);
    Matrix SingleSpatial(const Matrix& Aao, Matrix& Cmo);
}

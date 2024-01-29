#pragma once

#include "eigen.h"

namespace Numpy {
    Matrix<> Repeat(const Matrix<>& A, int count, int axis);
}

#pragma once

#include "integral.h"

namespace Population {
    Vector<> Mulliken(const System& system, const Integrals& ints, const Matrix<>& D);
}

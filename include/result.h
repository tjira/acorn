#pragma once

#include "eigen.h"

struct Result {
    struct RestrictedHartreeFock {
        double E; Vector<> eps; Matrix<> C, D, G, H;
    } rhf;
    struct RestrictedMollerPlesset {
        double Ecorr; Matrix<> G, H;
    } rmp;
    double Etot; Matrix<> G, H;
};

#pragma once

#include "eigen.h"

struct Result {
    struct RestrictedHartreeFock {
        double E; Vector<> eps; Matrix<> C, D, G, H;
    } rhf;
    struct RestrictedMollerPlesset {
        double Ecorr; Matrix<> G, H;
    } rmp;
    struct ModelSolver {
        std::vector<std::vector<Vector<std::complex<double>>>> states;
        Vector<> E, U, r, k;
    } msv;
    double Etot; Matrix<> G, H;
};

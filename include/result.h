#pragma once

#include "eigen.h"

struct Result {
    struct RestrictedHartreeFock {
        double E; Vector<> eps; Matrix<> C, D, G, H;
    } rhf;
    struct RestrictedConfigurationInteraction {
        Vector<> eps; Matrix<> C, F, G, H;
    } rci;
    struct RestrictedMollerPlesset {
        double Ecorr; Matrix<> G, H;
    } rmp;
    struct ModelSolver {
        std::vector<Vector<std::complex<double>>> optstates; 
        Vector<> opten; Vector<> r, k; Matrix<> U;
    } msv;
    double Etot; Matrix<> G, H;
};

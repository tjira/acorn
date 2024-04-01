#pragma once

#include "eigen.h"

struct Result {
    struct RestrictedHartreeFock {
        double E; Vector<> eps; Matrix<> C, D, G, H;
    } rhf;
    struct RestrictedConfigurationInteraction {
        Vector<> Eexc; Matrix<> C, F, G, H;
    } rci;
    struct RestrictedMollerPlesset {
        double Ecorr; Matrix<> G, H;
    } rmp;
    struct UnrestrictedHartreeFock {
        double E; Vector<> epsa, epsb; Matrix<> Ca, Cb, Da, Db, G, H;
    } uhf;
    struct UnrestrictedMollerPlesset {
        double Ecorr; Matrix<> G, H;
    } ump;
    struct ModelSolver {
        std::vector<Matrix<std::complex<double>>> optstates;
        Vector<> opten; Vector<> r, k, t, f; Matrix<> U;
        std::vector<Vector<std::complex<double>>> acfs;
        std::vector<Vector<double>> spectra;
    } msv;
    double Etot; Vector <> Eexc; Matrix<> G, H;
};

#include "restrictedhartreefock.h"

double RestrictedHartreeFock::energy(const System& system) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    libint2::initialize(); ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system); libint2::finalize();

    // run the restricted Hartree-Fock method and return the total energy
    return RestrictedHartreeFock(opt).run(system, ints, {}, false).Etot;
}

Method::Result RestrictedHartreeFock::run(const System& system, const Integrals& ints, Result, bool print) const {
    // create the antisymetrized Coulomb integral and define contraction axes with DIIS
    libint2::DIIS<Matrix<>> diis(2, 5); Eigen::IndexPair<int> first(2, 0), second(3, 1);
    Tensor<> ERI = ints.J - 0.5 * ints.J.shuffle(Eigen::array<int, 4>{0, 3, 2, 1});

    // define result struct and initialize the density matrix
    Result res; res.rhf.D = Matrix<>::Zero(ints.S.rows(), ints.S.cols());

    // define the Hamiltonian with Fock matrix and calculate energy
    Matrix<> H = ints.T + ints.V, F = H; res.rhf.E = 0.5 * res.rhf.D.cwiseProduct(H + F).sum();

    // print the iteration header
    if (print) std::printf("\nITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    for (int i = 1; i <= opt.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // calculate the Fock matrix
        F = H + toMatrix(ERI.contract(toTensor(res.rhf.D), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));

        // exrapolate the fock matrix
        Matrix<> e = ints.S * res.rhf.D * F - F * res.rhf.D * ints.S; diis.extrapolate(F, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix<>> solver(F, ints.S);

        // exteract the eigenvalues and eigenvectors
        res.rhf.C = solver.eigenvectors(), res.rhf.eps = solver.eigenvalues();

        // save previous density matrix and energy
        Matrix<> Dp = res.rhf.D; double Ep = res.rhf.E;

        // calculate the new density matrix
        res.rhf.D = 2 * res.rhf.C.leftCols(system.nocc()) * res.rhf.C.leftCols(system.nocc()).transpose();

        // calculate the new energy
        res.rhf.E = 0.5 * res.rhf.D.cwiseProduct(H + F).sum();

        // calculate the E and D errors
        double Eerr = std::abs(res.rhf.E - Ep), Derr = (res.rhf.D - Dp).norm();

        // extract the elapsed time
        std::string elapsed = Timer::Format(Timer::Elapsed(start));

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s\n", i, res.rhf.E, Eerr, Derr, elapsed.c_str());

        // finish if covergence reached
        if (Eerr < opt.thresh && Derr < opt.thresh) break;
        else if (i == opt.maxiter) {
            throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN SCF REACHED.");
        }
    }

    // assign total energy and return the struct
    res.Etot = res.rhf.E + system.repulsion(); return res;
}

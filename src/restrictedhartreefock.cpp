#include "restrictedhartreefock.h"

// include DIIS from libint
#include <libint2/diis.h>

std::tuple<Result, Integrals> RestrictedHartreeFock::run(const System& system, Result res, bool print) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system);

    // run the restricted Hartree-Fock method and return
    return {run(system, ints, res, print), ints};
}

Result RestrictedHartreeFock::run(const System& system, const Integrals& ints, Result res, bool print) const {
    // create the antisymmetrized Coulomb integral and define contraction axes with DIIS
    libint2::DIIS<Matrix<>> diis(1, 5); Eigen::IndexPair<int> first(2, 0), second(3, 1);
    Tensor<> ERI = ints.J - 0.5 * ints.J.shuffle(Eigen::array<int, 4>{0, 3, 2, 1});

    // initialize the density matrix
    res.rhf.D = res.rhf.D.size() ? res.rhf.D : Matrix<>::Zero(ints.S.rows(), ints.S.cols());

    // define the Hamiltonian and Fock matrix
    Matrix<> H = ints.T + ints.V, F = H;

    // print the iteration header
    if (print) std::printf("ITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    for (int i = 1; i <= opt.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // calculate the Fock matrix
        F = H + toMatrix(ERI.contract(toTensor(res.rhf.D), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));

        // exrapolate the fock matrix
        Matrix<> e = ints.S * res.rhf.D * F - F * res.rhf.D * ints.S; if (i > 1) diis.extrapolate(F, e);

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
        if (Eerr < opt.thresh && Derr < opt.thresh) {if (print) std::cout << std::endl; break;}
        else if (i == opt.maxiter) {
            throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN SCF REACHED.");
        }
    }

    // assign total energy and return the struct
    res.Etot = res.rhf.E + system.repulsion(); return res;
}

Result RestrictedHartreeFock::gradient(const System& system, const Integrals& ints, Result res, bool) const {
    // extract the useful stuff from the calculated integrals, define all the contraction axes and create the gradient matrix
    Tensor<3> dS1 = ints.dS.slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>({0, 0, 0}, {ints.dS.dimension(0), ints.dS.dimension(1), 3});
    Tensor<3> dT1 = ints.dT.slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>({0, 0, 0}, {ints.dT.dimension(0), ints.dT.dimension(1), 3});
    Tensor<3> dV1 = ints.dV.slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>({0, 0, 0}, {ints.dV.dimension(0), ints.dV.dimension(1), 3});
    Eigen::IndexPair<int> first(2, 0), second(3, 1), third(0, 0), fourth(1, 1); int nocc = system.nocc(); res.rhf.G = Matrix<>(system.getAtoms().size(), 3);

    // get the atom to shell map
    auto atom2shell = system.getShells().atom2shell(system.getAtoms<libint2::Atom>());

    // define the energy weighted density matrix
    Tensor<2> W(res.rhf.C.rows(), res.rhf.C.cols());

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * res.rhf.C.leftCols(nocc).row(i).cwiseProduct(res.rhf.C.leftCols(nocc).row(j)) * res.rhf.eps.topRows(nocc);
        }
    }

    // calculate the derivative of the ERI
    Tensor<3> dERI = (ints.dJ - 0.5 * ints.dJ.shuffle(Eigen::array<int, 5>{0, 3, 2, 1, 4})).contract(toTensor(res.rhf.D), Eigen::array<Eigen::IndexPair<int>, 2>{first, second});

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < res.rhf.G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += system.getShells().at(shell).size();

        // define the core Hamiltonian derivative, atomic slices for overlap tensor and density matrix
        Eigen::array<Eigen::Index, 3> Soff = {si, 0, 0}, Sext = {ss, res.rhf.D.cols(), 3}; Eigen::array<Eigen::Index, 2> Doff = {si, 0}, Dext = {ss, res.rhf.D.cols()};
        Tensor<3> dHcore = ints.dV.slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>({0, 0, 6 + i * 3}, {res.rhf.D.rows(), res.rhf.D.cols(), 3});

        // get the slice to add to the core Hamiltonian derivative
        auto HS = (dT1 + dV1).slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>(Soff, Sext);

        // add the slices to the core Hamiltonian derivative
        dHcore.slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>({0, si, 0}, {res.rhf.D.rows(), ss, 3}) += HS.shuffle(Eigen::array<Eigen::Index, 3>{1, 0, 2});
        dHcore.slice<Eigen::array<Eigen::Index, 3>, Eigen::array<Eigen::Index, 3>>({si, 0, 0}, {ss, res.rhf.D.cols(), 3}) += HS.shuffle(Eigen::array<Eigen::Index, 3>{0, 1, 2});

        // contract the tensors and add them to the gradient
        res.rhf.G.row(i) += 2 * toVector(dERI.slice(Soff, Sext).contract(toTensor(res.rhf.D).slice(Doff, Dext), Eigen::array<Eigen::IndexPair<int>, 2>{third, fourth}));
        res.rhf.G.row(i) -= 2 * toVector(dS1.slice(Soff, Sext).contract(W.slice(Doff, Dext), Eigen::array<Eigen::IndexPair<int>, 2>{third, fourth}));
        res.rhf.G.row(i) += toVector(dHcore.contract(toTensor(res.rhf.D), Eigen::array<Eigen::IndexPair<int>, 2>{third, fourth}));
    }

    // add nuclear contribution and return the gradient
    res.rhf.G += system.drepulsion(), res.G = res.rhf.G; return res;
}

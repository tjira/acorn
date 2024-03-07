#include "unrestrictedhartreefock.h"

// include DIIS from libint
#include <libint2/diis.h>

std::tuple<Result, Integrals> UnrestrictedHartreeFock::run(const System& system, Result res, bool print) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system);

    // run the unrestricted Hartree-Fock method and return
    return {run(system, ints, res, print), ints};
}

Result UnrestrictedHartreeFock::run(const System& system, const Integrals& ints, Result res, bool print) const {
    // create the antisymmetrized Coulomb integral and define contraction axes
    Tensor<> ERI = ints.J - 0.5 * ints.J.shuffle(Eigen::array<int, 4>{0, 3, 2, 1});

    // initialize DIIS and contraction axes
    libint2::DIIS<Matrix<>> diisa(1, 5), diisb(1, 5);
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // calculate number of alpha and beta electrons
    int na = (system.getElectrons() + system.getMulti() - 1) / 2;
    int nb = (system.getElectrons() - system.getMulti() + 1) / 2;

    // initialize the density matrices
    res.uhf.Da = res.uhf.Da.size() ? res.uhf.Da : Matrix<>::Zero(ints.S.rows(), ints.S.cols());
    res.uhf.Db = res.uhf.Db.size() ? res.uhf.Db : Matrix<>::Zero(ints.S.rows(), ints.S.cols());

    // calculate the exchange tensor
    Tensor<> K = ints.J.shuffle(Eigen::array<int, 4>{0, 3, 2, 1}); 

    // define the Hamiltonian and Fock matrices
    Matrix<> H = ints.T + ints.V, Fa = H, Fb = H;

    // print the iteration header
    if (print) std::printf("ITER       Eel [Eh]         |dE|    |dD_a|   |dD_b|      TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= opt.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // create coulomb and exchange matrices for both spins
        Matrix<> Ja = toMatrix(ints.J.contract(toTensor(res.uhf.Da), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));
        Matrix<> Jb = toMatrix(ints.J.contract(toTensor(res.uhf.Db), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));
        Matrix<> Ka = toMatrix(K.contract(toTensor(res.uhf.Da), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));
        Matrix<> Kb = toMatrix(K.contract(toTensor(res.uhf.Db), Eigen::array<Eigen::IndexPair<int>, 2>{first, second}));
        
        // calculate the Fock matrices
        Fa = H + 0.5 * (Ja + Jb - Ka), Fb = H + 0.5 * (Ja + Jb - Kb);

        // calculate the error vector
        Matrix<> e = ints.S * res.uhf.Da * Fa - Fa * res.uhf.Da * ints.S + ints.S * res.uhf.Db * Fb - Fb * res.uhf.Db * ints.S;

        // extrapolate the Fock matrices
        if (i > 1) diisa.extrapolate(Fa, e), diisb.extrapolate(Fb, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix<>> solvera(Fa, ints.S);
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix<>> solverb(Fb, ints.S);

        // exteract the eigenvalues and eigenvectors
        res.uhf.Ca = solvera.eigenvectors(), res.uhf.epsa = solvera.eigenvalues();
        res.uhf.Cb = solverb.eigenvectors(), res.uhf.epsb = solverb.eigenvalues();

        // save previous values of D and E
        Matrix<> Dap = res.uhf.Da, Dbp = res.uhf.Db; double Ep = res.uhf.E;

        // calculate the density matrices
        res.uhf.Da = 2 * res.uhf.Ca.leftCols(na) * res.uhf.Ca.leftCols(na).transpose();
        res.uhf.Db = 2 * res.uhf.Cb.leftCols(nb) * res.uhf.Cb.leftCols(nb).transpose();

        // calculate the new energy
        res.uhf.E = 0.25 * (res.uhf.Da.cwiseProduct(H + Fa).sum() + res.uhf.Db.cwiseProduct(H + Fb).sum());

        // calculate the E and D errors
        double Eerr = std::abs(res.uhf.E - Ep), Daerr = (res.uhf.Da - Dap).norm(), Dberr = (res.uhf.Db - Dbp).norm();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %.2e %s\n", i, res.uhf.E, Eerr, Daerr, Dberr, Timer::Format(Timer::Elapsed(start)).c_str());

        // finish if covergence reached
        if (Eerr < opt.thresh && Daerr < opt.thresh && Dberr < opt.thresh) {if (print) std::cout << std::endl; break;}
        else if (i == opt.maxiter) {
            throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN SCF REACHED.");
        }
    }

    // assign total energy and return the struct
    res.Etot = res.uhf.E + system.repulsion(); return res;
}

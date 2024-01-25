#include "modelsolver.h"
#include "modelsystem.h"

#include <iostream>

#define I std::complex<double>(0, 1)

Result ModelSolver::run(const ModelSystem& system, Result res, bool print) {
    if (system.potential.size() == 1) return runad(system, res, print);
    else return runnad(system, res, print);
}

Result ModelSolver::runad(const ModelSystem& system, Result res, bool print) {
    // define the real and Fourier space
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid);

    // fill the real space and calculate dx
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0).at(0) + (system.limits.at(0).at(1) - system.limits.at(0).at(0)) * i / (system.ngrid - 1);
    } double dx = res.msv.r(1) - res.msv.r(0);

    // fill the Fourier space
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / dx); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());

    // calculate the potential function values
    res.msv.U = Expression(system.potential.at(0).at(0)).eval(res.msv.r);

    // define the initial wavefunction
    Vector<std::complex<double>> psi = (-(res.msv.r.array() - 0.8).pow(2)).exp();

    // define the vector of optimal wavefunction and energies if imaginary time is selected
    if (!opt.real) res.msv.optstates = std::vector<Vector<std::complex<double>>>(opt.nstate, psi), res.msv.opten = Vector<>(opt.nstate); 

    // define the real and imaginary space operators
    Vector<std::complex<double>> K = (-0.5 * res.msv.k.array().pow(2) * opt.step).exp();
    Vector<std::complex<double>> R = (-0.5 * res.msv.U.array() * opt.step).exp();

    // real time dynamics operators
    if (opt.real) {
        K = (-0.5 * I * res.msv.k.array().pow(2) * opt.step).exp();
        R = (-0.5 * I * res.msv.U.array() * opt.step).exp();

        if (opt.optimize) {
            // create imaginary options and set some necessary values
            Options imopt = opt; imopt.real = false, imopt.iters = 100000;

            // optimize the wavefunctions
            res.msv.optstates = ModelSolver(imopt).run(system, res, false).msv.optstates;
        }
    }

    // print the newline
    if (print) std::printf("\n");

    // loop over all states
    for (int i = 0; i < opt.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n ITER        Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opt.real ? "(RT)" : "(IT)");

        // assign the guess to psi, create a vector that will hold the state dynamics and define vector of energies
        if (res.msv.optstates.size()) {psi = res.msv.optstates.at(i);} std::vector<Vector<std::complex<double>>> wfns = {psi};

        // calculate the total energy of the guess function
        Vector<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * EigenFourier(res.msv.k.array().pow(2) * EigenFourier(psi).array(), true).array();
        Vector<std::complex<double>> Ep = EigenConj(psi).array() * res.msv.U.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

        // vector of state energy evolution
        std::vector<double> energies = {E};

        // propagate the current state
        for (int j = 1; j <= opt.iters; j++) {
            // save the previous values
            Vector<std::complex<double>> Dprev = psi.array().abs2(); double Eprev = E;

            // Trotter formula
            psi = R.array() * psi.array();
            psi = EigenFourier(K.array() * EigenFourier(psi).array(), true);
            psi = R.array() * psi.array();

            // subtract lower eigenstates
            for (int k = 0; k < i && !opt.real; k++) {
                psi = psi.array() - (EigenConj(res.msv.optstates.at(k)).array() * psi.array()).sum() * dx * res.msv.optstates.at(k).array();
            }

            // normalize the wavefunction
            psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);

            // calculate the total energy
            Ek = 0.5 * EigenConj(psi).array() * EigenFourier(res.msv.k.array().pow(2) * EigenFourier(psi).array(), true).array();
            Ep = EigenConj(psi).array() * res.msv.U.array() * psi.array(); E = (Ek + Ep).sum().real() * dx;

            // calculate the errors
            double Eerr = std::abs(E - Eprev), Derr = (psi.array().abs2() - Dprev.array()).abs2().sum();

            // print the iteration
            if (print) std::printf("%6d %20.14f %.2e %.2e\n", j, E, Eerr, Derr);

            // push back the wfns and energies
            wfns.push_back(psi), energies.push_back(E);

            // end the loop if converged
            if (!opt.real && Eerr < opt.thresh && Derr < opt.thresh) break;
            else if (j == opt.iters && !opt.real) {
                throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN ITP REACHED.");
            }
        }

        // assign the energy and wavefunction
        if (!opt.real) res.msv.optstates.at(i) = psi, res.msv.opten(i) = E;

        // save the state wavefunction
        ModelSystem::SaveWavefunction("state" + std::to_string(i) + ".dat", res.msv.r, wfns, energies);
    }

    // return the results
    return res;
}

Result ModelSolver::runnad(const ModelSystem& system, Result res, bool print) {
    // define the real space
    res.msv.r = Vector<>(system.ngrid), res.msv.U = Matrix<>(system.ngrid, system.potential.size());

    // energies container, wfn vector and wfn container and transformation matrix container
    std::vector<double> energies; Vector<std::complex<double>> psi(2 * system.ngrid);
    std::vector<Matrix<std::complex<double>>> chis; std::vector<Matrix<>> transform;

    // fill the real space and calculate dx
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0).at(0) + (system.limits.at(0).at(1) - system.limits.at(0).at(0)) * i / (system.ngrid - 1);
    } double dx = res.msv.r(1) - res.msv.r(0);

    // define the adiabatic kinetic matrix and Hamiltonian
    Matrix<> T(system.ngrid, system.ngrid), H(2 * system.ngrid, 2 * system.ngrid);

    // fill the kinetic matrix
    for (int i = 0; i < system.ngrid; i++) {
        for (int j = 0; j < system.ngrid; j++) {
            T(i, j) = 1.0 / (2 * system.m * dx * dx) * std::pow(-1, i + j);
            T(i, j) *= i == j ? M_PI * M_PI / 3 : 2.0 / (i - j) / (i - j);
        }
    }

    // fill the Hamiltonian with the kinetic matrix
    H.block(0, 0, system.ngrid, system.ngrid) = T, H.block(system.ngrid, system.ngrid, system.ngrid, system.ngrid) = T;

    // create the potential energy expressions
    std::vector<std::vector<Expression>> uexpr = {
        {Expression(system.potential.at(0).at(0)), Expression(system.potential.at(0).at(1))},
        {Expression(system.potential.at(1).at(0)), Expression(system.potential.at(1).at(1))}
    };

    for (int i = 0; i < system.ngrid; i++) {
        // define the potential matrix
        Matrix<> V(system.potential.size(), system.potential.size());

        // fill the potential matrix
        for (int j = 0; j < V.rows(); j++) {
            for (int k = 0; k < V.cols(); k++) {
                V(j, k) = uexpr.at(j).at(k).eval(res.msv.r(i));
            }
        }

        // add the potential energy to the Hamiltonian
        H(i, i) += V(0, 0), H(i + system.ngrid, i + system.ngrid) += V(1, 1);
        H(i, i + system.ngrid) += V(0, 1), H(i + system.ngrid, i) += V(1, 0);

        // find the adiabatic basis
        Eigen::SelfAdjointEigenSolver<Matrix<>> solver(V); res.msv.U.row(i) = solver.eigenvalues();

        // push back the transformation matrix
        transform.push_back(solver.eigenvectors());
    }

    // create the hamiltonian solver in DVR basis
    Eigen::SelfAdjointEigenSolver<Matrix<>> solver(H);

    // extract the eigenvalues and eigenvectors
    Vector<> eigv = solver.eigenvalues(); Matrix<> eigf = solver.eigenvectors();

    // fill in the initial wavefunction in diabatic basis
    psi.head(system.ngrid) = (-(res.msv.r.array() + 15).pow(2)).exp() * (I * std::sqrt(0.06 * system.m) * res.msv.r.array()).exp();

    // normalize the wavefunction
    psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);

    // calculate DVR expansion coeficients
    Vector<std::complex<double>> ci(psi.size());
    for (int j = 0; j < psi.size(); j++) {
        ci(j) = (psi.transpose() * eigf.col(j))(0, 0) * dx;
    }

    // calculate the total energy
    double E = (EigenConj(ci).transpose() * H * ci)(0, 0).real();


    // print the iteration header
    if (print) std::printf("\n ITER        Eel [Eh]         |dE|     |dD|\n");

    for (int i = 0; i < opt.iters; i++) {
        // define the next wavefunction
        Vector<std::complex<double>> nextpsi(psi.size());

        // propagate the wavefunction
        for (int j = 0; j < psi.size(); j++) {
            nextpsi += ci(j) * eigf.col(j) * std::exp(-I * eigv(j) * std::complex<double>(i, 0) * opt.step);
        } psi = nextpsi;

        // normalize the wavefunction
        psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);

        // prepare the matrices for the adiabatic transform
        Matrix<std::complex<double>> psimat(system.ngrid, system.potential.size());
        Matrix<std::complex<double>> chi(system.ngrid, system.potential.size());

        // fill the matrix of wavefunction states as columns
        psimat << psi.head(system.ngrid), psi.tail(system.ngrid);

        // transform the wavefunction to adiabatic basis
        for (int j = 0; j < system.ngrid; j++) {
            chi.row(j) = transform.at(j).transpose() * psimat.row(j).transpose();
        }

        // save the wavefunction to the vector
        chis.push_back(chi); energies.push_back(E);

        // print the iteration
        if (print) std::printf("%6d %20.14f %.2e %.2e\n", i + 1, E, 0.0, 0.0);
    }

    // save the wfn dynamics
    for (int i = 0; i < system.potential.size(); i++) {
        std::vector<Vector<std::complex<double>>> state; for (auto chi : chis) state.push_back(chi.col(i));
        ModelSystem::SaveWavefunction("state" + std::to_string(i) + ".dat", res.msv.r, state, energies);
    }

    // return
    return res;
}

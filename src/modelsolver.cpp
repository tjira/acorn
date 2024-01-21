#include "modelsolver.h"

#define I std::complex<double>(0, 1)

Result ModelSolver::run(const ModelSystem& system, Result res, bool print) {
    // define the real and Fourier space
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid);

    // fill the real space and calculate dx
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0).at(0) + (system.limits.at(0).at(1) - system.limits.at(0).at(0)) * i / (system.ngrid - 1);
    } double dx = res.msv.r(1) - res.msv.r(0);

    // fill the Fourier space
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / dx); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());

    // calculate the potential function values
    res.msv.U = Expression(system.potential).eval(res.msv.r);

    // define the initial wavefunction and the energy vector
    Vector<std::complex<double>> psi = (-(res.msv.r.array() - 0.8).pow(2)).exp(); res.msv.E = Vector<>(opt.nstate);

    // define the vector of all wavefunctions
    res.msv.states = std::vector<std::vector<Vector<std::complex<double>>>>(opt.nstate, {psi}); 

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
            res.msv.states = ModelSolver(imopt).run(system, res, false).msv.states;

            // delete optimization steps
            for (size_t i = 0; i < res.msv.states.size(); i++) {
                res.msv.states.at(i) = {res.msv.states.at(i).at(res.msv.states.at(i).size() - 1)};
            }
        }
    }

    // calculate the total energy
    Vector<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * EigenFourier(res.msv.k.array().pow(2) * EigenFourier(psi).array(), true).array();
    Vector<std::complex<double>> Ep = EigenConj(psi).array() * res.msv.U.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

    // print the header
    if (print) std::printf("\n");

    // loop over all states
    for (int i = 0; i < opt.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n ITER        Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opt.real ? "(RT)" : "(IT)");

        // assign the psi WFN to the correct state
        psi = res.msv.states.at(i).at(res.msv.states.at(i).size() - 1);

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
                psi = psi.array() - (EigenConj(res.msv.states.at(k).at(res.msv.states.at(k).size() - 1)).array() * psi.array()).sum() * dx * res.msv.states.at(k).at(res.msv.states.at(k).size() - 1).array();
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

            // append the wave function
            res.msv.states.at(i).push_back(psi);

            // end the loop if converged
            if (!opt.real && Eerr < opt.thresh && Derr < opt.thresh) break;
            else if (j == opt.iters && !opt.real) {
                throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN ITP REACHED.");
            }
        }

        // assign the energy
        res.msv.E(i) = E;
    }

    // print the new line and return the results
    if (print) {std::printf("\n");} return res;
}

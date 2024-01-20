#include "modelsolver.h"

#define I std::complex<double>(0, 1)

ModelSolver::Result ModelSolver::run(const ModelSystem& system, Result res, bool print) {
    // define the real and Fourier space
    Vector<std::complex<double>> x(system.ngrid), k(system.ngrid);

    // fill the real space and calculate dx
    for (int i = 0; i < system.ngrid; i++) {
        x(i) = system.limits.at(0).at(0) + (system.limits.at(0).at(1) - system.limits.at(0).at(0)) * i / (system.ngrid - 1);
    } double dx = (x(1) - x(0)).real();

    // fill the Fourier space
    k.fill(2 * M_PI / k.size() / dx); for (int i = 0; i < k.size(); i++) k(i) *= i - (i < x.size() / 2 ? 0 : x.size());

    // calculate the potential function values
    Vector<std::complex<double>> V = Expression(system.potential).eval(x.real());

    // define the initial wavefunction and the energy vector
    Vector<std::complex<double>> psi = (-(x.array() - 0.8).pow(2)).exp(); Vector<> energy(opt.nstate);

    // define the vector of all wavefunctions
    std::vector<std::vector<Vector<std::complex<double>>>> states(opt.nstate, {psi}); 

    // define the real and imaginary space operators
    Vector<std::complex<double>> K = (-0.5 * k.array().pow(2) * opt.step).exp();
    Vector<std::complex<double>> R = (-0.5 * V.array() * opt.step).exp();

    // real time dynamics operators
    if (opt.real) {
        K = (-0.5 * I * k.array().pow(2) * opt.step).exp();
        R = (-0.5 * I * V.array() * opt.step).exp();

        if (opt.optimize) {
            // create imaginary options and set some necessary values
            Options imopt = opt; imopt.real = false, imopt.iters = 100000;

            // optimize the wavefunctions
            states = ModelSolver(imopt).run(system, res, false).states;

            // delete optimization steps
            for (size_t i = 0; i < states.size(); i++) {
                states.at(i) = {states.at(i).at(states.at(i).size() - 1)};
            }
        }
    }

    // calculate the total energy
    Vector<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * EigenFourier(k.array().pow(2) * EigenFourier(psi).array(), true).array();
    Vector<std::complex<double>> Ep = EigenConj(psi).array() * V.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

    // print the header
    if (print) std::printf("\n");

    // loop over all states
    for (int i = 0; i < opt.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n ITER        Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opt.real ? "(RT)" : "(IT)");

        // assign the psi WFN to the correct state
        psi = states.at(i).at(states.at(i).size() - 1);

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
                psi = psi.array() - (EigenConj(states.at(k).at(states.at(k).size() - 1)).array() * psi.array()).sum() * dx * states.at(k).at(states.at(k).size() - 1).array();
            }

            // normalize the wavefunction
            psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);

            // calculate the total energy
            Ek = 0.5 * EigenConj(psi).array() * EigenFourier(k.array().pow(2) * EigenFourier(psi).array(), true).array();
            Ep = EigenConj(psi).array() * V.array() * psi.array(); E = (Ek + Ep).sum().real() * dx;

            // calculate the errors
            double Eerr = std::abs(E - Eprev), Derr = (psi.array().abs2() - Dprev.array()).abs2().sum();

            // print the iteration
            if (print) std::printf("%6d %20.14f %.2e %.2e\n", j, E, Eerr, Derr);

            // append the wave function
            states.at(i).push_back(psi);

            // end the loop if converged
            if (!opt.real && Eerr < opt.thresh && Derr < opt.thresh) break;
            else if (j == opt.iters && !opt.real) {
                throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN ITP REACHED.");
            }
        }

        // assign the energy
        energy(i) = E;
    }

    // print the new line and return the results
    if (print) {std::printf("\n");} return {states, energy, x.real(), V.real()};
}

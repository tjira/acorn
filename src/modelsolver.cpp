#include "modelsolver.h"
#include "modelsystem.h"

#define I std::complex<double>(0, 1)

Result ModelSolver::run(const ModelSystem& system, Result res, bool print) {
    if (system.potential.size() == 1) return runad(system, res, print);
    else return runnad(system, res, print);
}

Result ModelSolver::runad(const ModelSystem& system, Result res, bool print) {
    // define the real and Fourier space along with time and frequency domains
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid), res.msv.t = Vector<>(opta.iters + 1), res.msv.f = Vector<>(opta.iters + 1);

    // fill the real space and calculate dx
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } double dx = res.msv.r(1) - res.msv.r(0);

    // fill the time domain
    for (int i = 0; i <= opta.iters; i++) res.msv.t(i) = i * opta.step;

    // fill the frequency and momentum domain
    res.msv.f.fill(2 * M_PI / res.msv.f.size() / opta.step); for (int i = 0; i < res.msv.f.size(); i++) res.msv.f(i) *= i - (i < res.msv.f.size() / 2 ? 0 : res.msv.f.size());
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / dx); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());

    // calculate the potential function values
    res.msv.U = Expression(system.potential.at(0).at(0)).eval(res.msv.r);

    // define the initial wavefunction and normalize it
    Vector<std::complex<double>> psi = Expression(opta.guess).eval(res.msv.r); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);

    // define the vector of optimal wavefunction and energies if imaginary time is selected
    if (!opta.real) res.msv.optstates = std::vector<Vector<std::complex<double>>>(opta.nstate, psi), res.msv.opten = Vector<>(opta.nstate); 

    // define the real and imaginary space operators
    Vector<std::complex<double>> K = (-0.5 * res.msv.k.array().pow(2) * opta.step / opta.mass).exp();
    Vector<std::complex<double>> R = (-0.5 * res.msv.U.array() * opta.step).exp();

    // real time dynamics operators
    if (opta.real) {
        K = (-0.5 * I * res.msv.k.array().pow(2) * opta.step / opta.mass).exp();
        R = (-0.5 * I * res.msv.U.array() * opta.step).exp();

        if (opta.optimize) {
            // create imaginary options and set some necessary values
            OptionsAdiabatic imopt = opta; imopt.real = false, imopt.iters = 100000;

            // optimize the wavefunctions
            res.msv.optstates = ModelSolver(imopt).run(system, res, false).msv.optstates;
        }
    }

    // loop over all states
    for (int i = 0; i < opta.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n ITER        Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opta.real ? "(RT)" : "(IT)");

        // assign the guess to psi, create a vector that will hold the state dynamics and define vector of energies
        if (res.msv.optstates.size()) {psi = res.msv.optstates.at(i);} std::vector<Vector<std::complex<double>>> wfns = {psi};

        // calculate the total energy of the guess function
        Vector<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * EigenFourier(res.msv.k.array().pow(2) * EigenFourier(psi).array(), true).array();
        Vector<std::complex<double>> Ep = EigenConj(psi).array() * res.msv.U.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

        // vector of state energy evolution
        std::vector<double> energies = {E};

        // propagate the current state
        for (int j = 1; j <= opta.iters; j++) {
            // save the previous values
            Vector<std::complex<double>> Dprev = psi.array().abs2(); double Eprev = E;

            // Trotter formula
            psi = R.array() * psi.array();
            psi = EigenFourier(K.array() * EigenFourier(psi).array(), true);
            psi = R.array() * psi.array();

            // subtract lower eigenstates
            for (int k = 0; k < i && !opta.real; k++) {
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
            if (!opta.real && Eerr < opta.thresh && Derr < opta.thresh) break;
            else if (j == opta.iters && !opta.real) {
                throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN ITP REACHED.");
            }
        }

        // assign the energy and wavefunction
        if (!opta.real) res.msv.optstates.at(i) = psi, res.msv.opten(i) = E;

        // block to calculate spectrum
        if (opta.real && opta.spectrum) {
            // initialize the autocorrelation function
            Vector<std::complex<double>> acf(wfns.size());

            // calculate the autocorrelation function
            for (size_t j = 0; j < wfns.size(); j++) {
                acf(j) = (wfns.at(0).array() * EigenConj(wfns.at(j)).array()).sum() * dx;
            }

            // append and perform the fourier transform
            res.msv.acfs.push_back(acf); res.msv.spectra.push_back(EigenFourier(acf));
        }

        // save the state wavefunction
        ModelSystem::SaveWavefunction(opta.folder + "/state" + std::to_string(i) + ".dat", res.msv.r, wfns, energies);
    }

    // print the newline
    if (print && !opta.real) std::printf("\n");

    // return the results
    return res;
}

Result ModelSolver::runnad(const ModelSystem& system, Result res, bool print) {
    // define the real and momentum space
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid);

    // energies container, wfn vector and wfn array
    std::vector<double> energies; Matrix<std::complex<double>> psi(system.ngrid, 2); std::vector<Matrix<std::complex<double>>> psis;

    // fill the real space and calculate dx
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } double dx = res.msv.r(1) - res.msv.r(0);

    // fill the momentum space
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / dx); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());

    // fill in the initial wavefunction
    for (size_t i = 0; i < optn.guess.size(); i++) {
        psi.block(0, i, system.ngrid, 1) = (Expression(optn.guess.at(i)).eval(res.msv.r).array() * (I * std::sqrt(0.06 * system.m) * res.msv.r.array()).exp());
    }

    // normalize the wavefunction
    psi.block(0, 0, system.ngrid, 1) = psi.block(0, 0, system.ngrid, 1).array() / std::sqrt(psi.block(0, 0, system.ngrid, 1).array().abs2().sum() * dx);

    // create the potential energy expressions
    std::vector<std::vector<Expression>> uexpr = {
        {Expression(system.potential.at(0).at(0)), Expression(system.potential.at(0).at(1))},
        {Expression(system.potential.at(1).at(0)), Expression(system.potential.at(1).at(1))}
    }; res.msv.U = Matrix<>(system.ngrid, system.potential.size());

    // fill the potential
    res.msv.U.col(0) = uexpr.at(0).at(0).eval(res.msv.r);
    res.msv.U.col(1) = uexpr.at(1).at(1).eval(res.msv.r);

    // append the first wfn and energy
    psis.push_back(psi), energies.push_back(0);

    // define the real and imaginary space operators
    std::vector<std::vector<Vector<std::complex<double>>>> K(2, std::vector<Vector<std::complex<double>>>(2)); auto R = K;

    // fill the kinetic propagator
    K.at(0).at(0) = (-0.5 * I * res.msv.k.array().pow(2) * optn.step / optn.mass).exp();
    K.at(1).at(1) = (-0.5 * I * res.msv.k.array().pow(2) * optn.step / optn.mass).exp();

    // fill the potential propagator
    Vector<std::complex<double>> D = 4 * uexpr.at(1).at(0).eval(res.msv.r).array().abs().pow(2) + (uexpr.at(0).at(0).eval(res.msv.r) - uexpr.at(1).at(1).eval(res.msv.r)).array().pow(2);
    Vector<std::complex<double>> a = (-0.25 * I * (uexpr.at(0).at(0).eval(res.msv.r) + uexpr.at(1).at(1).eval(res.msv.r)) * optn.step).array().exp();
    Vector<std::complex<double>> b = (0.25 * D.array().sqrt() * optn.step).cos();
    Vector<std::complex<double>> c = I * (0.25 * D.array().sqrt() * optn.step).sin() / D.array().sqrt();
    R.at(0).at(0) = a.array() * (b.array() + c.array() * (uexpr.at(1).at(1).eval(res.msv.r) - uexpr.at(0).at(0).eval(res.msv.r)).array());
    R.at(0).at(1) = -2 * a.array() * c.array() * uexpr.at(0).at(1).eval(res.msv.r).array();
    R.at(1).at(0) = -2 * a.array() * c.array() * uexpr.at(0).at(1).eval(res.msv.r).array();
    R.at(1).at(1) = a.array() * (b.array() + c.array() * (uexpr.at(0).at(0).eval(res.msv.r) - uexpr.at(1).at(1).eval(res.msv.r)).array());

    // print the iteration header
    if (print) std::printf(" ITER        Eel [Eh]         |dE|     |dD|\n");

    // propagate the wavefunction
    for (int i = 0; i < optn.iters; i++) {
        psi.col(0) = EigenFourier(R.at(0).at(0).array() * psi.col(0).array() + R.at(0).at(1).array() * psi.col(1).array());
        psi.col(1) = EigenFourier(R.at(1).at(0).array() * psi.col(0).array() + R.at(1).at(1).array() * psi.col(1).array());
        psi.col(0) = EigenFourier(K.at(0).at(0).array() * psi.col(0).array(), true);
        psi.col(1) = EigenFourier(K.at(1).at(1).array() * psi.col(1).array(), true);
        psi.col(0) = (R.at(0).at(0).array() * psi.col(0).array() + R.at(0).at(1).array() * psi.col(1).array());
        psi.col(1) = (R.at(1).at(0).array() * psi.col(0).array() + R.at(1).at(1).array() * psi.col(1).array());

        // calculate the errors
        double E = 0, Eerr = 0, Derr = (psi.array().abs2() - psis.at(psis.size() - 1).array().abs2()).abs2().sum();

        // append the wfn and energy
        psis.push_back(psi), energies.push_back(0);

        // print the iteration
        if (print) std::printf("%6d %20.14f %.2e %.2e\n", i + 1, E, Eerr, Derr);
    }

    // save the wfn dynamics
    for (size_t i = 0; i < system.potential.size(); i++) {
        std::vector<Vector<std::complex<double>>> state; for (auto psi : psis) state.push_back(psi.col(i));
        ModelSystem::SaveWavefunction(optn.folder + "/state" + std::to_string(i) + ".dat", res.msv.r, state, energies);
    }

    // return
    return res;
}

#include "modelsolver.h"
#include "modelsystem.h"

#define I std::complex<double>(0, 1)

Result ModelSolver::run(const ModelSystem& system, Result res, bool print) {
    if (system.potential.size() == 1) return runad(system, res, print);
    else return runnad(system, res, print);
}

Result ModelSolver::runad(const ModelSystem& system, Result res, bool print) {
    // obtain the simulation dimension
    int dim; dim = StringContains(system.potential.at(0).at(0), 'y') ? 2 : 1;

    // define the real and Fourier space along with time and frequency domains and wavefunction
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid), res.msv.t = Vector<>(opta.iters + 1), res.msv.f = Vector<>(opta.iters + 1); Matrix<std::complex<double>> psi;

    // fill the real space and calculate dr
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } double dr = std::pow(res.msv.r(1) - res.msv.r(0), dim);

    // fill the time domain
    for (int i = 0; i <= opta.iters; i++) res.msv.t(i) = i * opta.step;

    // fill the momentum and frequency domain
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / (res.msv.r(1) - res.msv.r(0))); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());
    res.msv.f.fill(2 * M_PI / res.msv.f.size() / opta.step); for (int i = 0; i < res.msv.f.size(); i++) res.msv.f(i) *= i - (i < res.msv.f.size() / 2 ? 0 : res.msv.f.size());

    // define the 2-dimensional variables and potential values
    Matrix<> x(system.ngrid, system.ngrid), y(system.ngrid, system.ngrid), k(system.ngrid, system.ngrid), l(system.ngrid, system.ngrid), U;

    // fill the 2-dimensional variables
    for (int i = 0; i < system.ngrid; i++) {
        y.col(i) = res.msv.r, l.col(i) = res.msv.k;
    } x = y.transpose(), k = l.transpose();
    
    // initialize the potential matrix
    res.msv.U = Matrix<>((int)std::pow(system.ngrid, dim), 1);

    // calculate the potential values
    if (dim == 1) U = Expression(system.potential.at(0).at(0)).eval(res.msv.r);
    if (dim == 2) U = Expression(system.potential.at(0).at(0)).eval(x, y);

    // define the initial wavefunction and normalize it
    if (dim == 1) {psi = Expression(opta.guess).eval(res.msv.r).col(0); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dr);}
    if (dim == 2) {psi = Expression(opta.guess).eval(x, y); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dr);}

    // assign correct variables for 1-dimensional model
    if (dim == 1) x = res.msv.r, k = res.msv.k, y = Matrix<>::Zero(system.ngrid, 1), l = Matrix<>::Zero(system.ngrid, 1);

    // fill the potential matrix
    for (int i = 0; i < res.msv.U.rows(); i++) {
        res.msv.U(i, 0) = U(i % system.ngrid, i / system.ngrid);
    }

    // define the vector of optimal wavefunction and energies if imaginary time is selected
    if (!opta.real) res.msv.optstates = std::vector<Matrix<std::complex<double>>>(opta.nstate, psi), res.msv.opten = Vector<>(opta.nstate); 

    // define the imaginary time operators
    Matrix<std::complex<double>> K = (-0.5 * (k.array().pow(2) + l.array().pow(2)) * opta.step / system.mass()).exp();
    Matrix<std::complex<double>> R = (-0.5 * U.array() * opta.step).exp();

    // real time dynamics operators
    if (opta.real) {
        if (opta.optimize) {
            // create imaginary options and with necessary values
            OptionsAdiabatic imopt = opta; imopt.real = false, imopt.iters = 1000;

            // optimize the wavefunctions
            auto optres = ModelSolver(imopt).run(system, res, false); res.msv.optstates = optres.msv.optstates, res.msv.opten = optres.msv.opten;
        }

        // change the potential for te real time dynamics
        if (dim == 1) U = Expression(opta.spectrum.potential).eval(res.msv.r);
        if (dim == 2) U = Expression(opta.spectrum.potential).eval(x, y);

        // define the real time operators
        K = (-0.5 * I * (k.array().pow(2) + l.array().pow(2)) * opta.step / system.mass()).exp();
        R = (-0.5 * I * U.array() * opta.step).exp();
    }

    // loop over all states
    for (int i = 0; i < opta.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n ITER        Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opta.real ? "(RT)" : "(IT)");

        // assign the guess to psi, create a vector that will hold the state dynamics and define vector of energies
        if (res.msv.optstates.size()) {psi = res.msv.optstates.at(i);} std::vector<Matrix<std::complex<double>>> wfns = {psi};

        // calculate the total energy of the guess function
        Matrix<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * Numpy::FFT((k.array().pow(2) + l.array().pow(2)) * Numpy::FFT(psi).array(), 1).array();
        Matrix<std::complex<double>> Ep = EigenConj(psi).array() * U.array() * psi.array(); double E = (Ek / system.mass() + Ep).sum().real() * dr;

        // vector of state energy evolution
        std::vector<double> energies = {E};

        // propagate the current state
        for (int j = 1; j <= opta.iters; j++) {
            // save the previous values
            Matrix<std::complex<double>> Dprev = psi.array().abs2(); double Eprev = E;

            // Trotter formula
            psi = R.array() * psi.array();
            psi = Numpy::FFT(K.array() * Numpy::FFT(psi).array(), 1);
            psi = R.array() * psi.array();

            // subtract lower eigenstates
            for (int k = 0; k < i && !opta.real; k++) {
                psi = psi.array() - (EigenConj(res.msv.optstates.at(k)).array() * psi.array()).sum() * dr * res.msv.optstates.at(k).array();
            }

            // normalize the wavefunction
            psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dr);

            // calculate the total energy
            Ek = 0.5 * EigenConj(psi).array() * Numpy::FFT((k.array().pow(2) + l.array().pow(2)) * Numpy::FFT(psi).array(), 1).array();
            Ep = EigenConj(psi).array() * U.array() * psi.array(); E = (Ek / system.mass() + Ep).sum().real() * dr;

            // calculate the errors
            double Eerr = std::abs(E - Eprev), Derr = (psi.array().abs2() - Dprev.array()).abs2().sum();

            // print the iteration
            if (print) std::printf("%6d %20.14f %.2e %.2e\n", j, E, Eerr, Derr);

            // push back the wfns and energies
            wfns.push_back(psi), energies.push_back(E);
        }

        // assign the energy and wavefunction
        if (!opta.real) res.msv.optstates.at(i) = psi, res.msv.opten(i) = E;

        // block to calculate spectrum
        if (opta.real && !opta.spectrum.potential.empty()) {
            // initialize the autocorrelation function
            Vector<std::complex<double>> acf(wfns.size());

            // calculate the autocorrelation function
            for (size_t j = 0; j < wfns.size(); j++) {
                acf(j) = (wfns.at(0).array() * EigenConj(wfns.at(j)).array()).sum() * dr;
            }

            // define the windowing function
            auto wfun = (0.5 - 0.5 * (2 * M_PI * res.msv.t.array() / (wfns.size() - 1)).cos());

            // append and perform the Fourier transform
            res.msv.acfs.push_back(acf); res.msv.spectra.push_back(Numpy::FFT(wfun.array() * acf.array()));

            // subtract the ground state energy to get the correct spectrum
            if (res.msv.opten.size()) {
                res.msv.f = res.msv.f.array() - res.msv.opten(0);
            }
        }

        // save the state wavefunction
        if (dim == 2) ModelSystem::SaveWavefunction(opta.folder + "/state" + std::to_string(i) + ".dat", x.row(0), y.col(0), wfns, energies);
        if (dim == 1) ModelSystem::SaveWavefunction(opta.folder + "/state" + std::to_string(i) + ".dat", x, wfns, energies);
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

    // fill the real space and calculate dr
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } double dr = res.msv.r(1) - res.msv.r(0);

    // fill the momentum space
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / dr); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());

    // fill in the initial wavefunction
    for (size_t i = 0; i < optn.guess.size(); i++) {
        psi.block(0, i, system.ngrid, 1) = (Expression(optn.guess.at(i)).eval(res.msv.r).array() * (I * std::sqrt(0.06 * system.m) * res.msv.r.array()).exp());
    }

    // normalize the wavefunction
    psi.block(0, 0, system.ngrid, 1) = psi.block(0, 0, system.ngrid, 1).array() / std::sqrt(psi.block(0, 0, system.ngrid, 1).array().abs2().sum() * dr);

    // create the potential
    std::vector<std::vector<Vector<std::complex<double>>>> V = {
        {Expression(system.potential.at(0).at(0)).eval(res.msv.r), Expression(system.potential.at(0).at(1)).eval(res.msv.r)},
        {Expression(system.potential.at(1).at(0)).eval(res.msv.r), Expression(system.potential.at(1).at(1)).eval(res.msv.r)}
    }; res.msv.U = Matrix<>(system.ngrid, system.potential.size());

    // save the potential to results
    res.msv.U.col(0) = Expression(system.potential.at(0).at(0)).eval(res.msv.r);
    res.msv.U.col(1) = Expression(system.potential.at(1).at(1)).eval(res.msv.r);

    // define the complex absorbing potential
    Vector<std::complex<double>> CAP = 0.1 * I * ((-0.1 * (res.msv.r.array() - system.limits.at(0)).pow(2)).exp() + (-0.1 * (res.msv.r.array() - system.limits.at(1)).pow(2)).exp());

    // subtract the complex absorbing potential
    V.at(0).at(0) = V.at(0).at(0).array() - CAP.array();
    V.at(1).at(1) = V.at(1).at(1).array() - CAP.array();

    // append the first wfn and energy
    psis.push_back(psi), energies.push_back(0);

    // define the real and imaginary space operators
    std::vector<std::vector<Vector<std::complex<double>>>> K(2, std::vector<Vector<std::complex<double>>>(2)); auto R = K;

    // fill the kinetic propagator
    K.at(0).at(0) = (-0.5 * I * res.msv.k.array().pow(2) * optn.step / system.mass()).exp();
    K.at(1).at(1) = (-0.5 * I * res.msv.k.array().pow(2) * optn.step / system.mass()).exp();

    // fill the potential propagator
    Vector<std::complex<double>> D = 4 * V.at(1).at(0).array().abs().pow(2) + (V.at(0).at(0) - res.msv.U.col(1)).array().pow(2);
    Vector<std::complex<double>> a = (-0.25 * I * (V.at(0).at(0) + res.msv.U.col(1)) * optn.step).array().exp();
    Vector<std::complex<double>> b = (0.25 * D.array().sqrt() * optn.step).cos();
    Vector<std::complex<double>> c = I * (0.25 * D.array().sqrt() * optn.step).sin() / D.array().sqrt();
    R.at(0).at(0) = a.array() * (b.array() + c.array() * (V.at(1).at(1) - V.at(0).at(0)).array());
    R.at(0).at(1) = -2 * a.array() * c.array() * V.at(0).at(1).array();
    R.at(1).at(0) = -2 * a.array() * c.array() * V.at(0).at(1).array();
    R.at(1).at(1) = a.array() * (b.array() + c.array() * (V.at(0).at(0) - V.at(1).at(1)).array());

    // print the iteration header
    if (print) std::printf(" ITER        Eel [Eh]         |dE|     |dD|\n");

    // propagate the wavefunction
    for (int i = 0; i < optn.iters; i++) {
        // apply the propagator
        psi.col(0) = Numpy::FFT(R.at(0).at(0).array() * psi.col(0).array() + R.at(0).at(1).array() * psi.col(1).array());
        psi.col(1) = Numpy::FFT(R.at(1).at(0).array() * psi.col(0).array() + R.at(1).at(1).array() * psi.col(1).array());
        psi.col(0) = Numpy::FFT(K.at(0).at(0).array() * psi.col(0).array(), 1);
        psi.col(1) = Numpy::FFT(K.at(1).at(1).array() * psi.col(1).array(), 1);
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
        std::vector<Matrix<std::complex<double>>> state; for (auto psi : psis) state.push_back(psi.col(i));
        ModelSystem::SaveWavefunction(optn.folder + "/state" + std::to_string(i) + ".dat", res.msv.r, state, energies);
    }

    // return
    return res;
}

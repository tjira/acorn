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

    // define the real and Fourier space along with time and frequency domains
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid), res.msv.t = Vector<>(opta.iters + 1), res.msv.f = Vector<>(opta.iters + opta.spectrum.zeropad + 1);

    // define the wavefunction, potential and independent variables container
    Matrix<std::complex<double>> psi; Matrix<> U, x, y, k, l;

    // fill the real space and calculate dr
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } double dr = std::pow(res.msv.r(1) - res.msv.r(0), dim);

    // fill the time domain
    for (int i = 0; i <= opta.iters; i++) res.msv.t(i) = i * opta.step;

    // fill the momentum and frequency domain
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / (res.msv.r(1) - res.msv.r(0))); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());
    res.msv.f.fill(2 * M_PI / res.msv.f.size() / (res.msv.t(1) - res.msv.t(0))); for (int i = 0; i < res.msv.f.size(); i++) res.msv.f(i) *= i - (i < res.msv.f.size() / 2 ? 0 : res.msv.f.size());

    // initialize the 2-dimensional variables
    if (dim == 2) {
        x = Matrix<>::Zero(system.ngrid, system.ngrid), y = Matrix<>::Zero(system.ngrid, system.ngrid);
        k = Matrix<>::Zero(system.ngrid, system.ngrid), l = Matrix<>::Zero(system.ngrid, system.ngrid);
    } else if (dim == 1) {
        x = Vector<>::Zero(system.ngrid), y = Vector<>::Zero(system.ngrid);
        k = Vector<>::Zero(system.ngrid), l = Vector<>::Zero(system.ngrid);
    }

    // fill the independent variables
    if (dim == 2) {
        for (int i = 0; i < system.ngrid; i++) {
            x.col(i) = res.msv.r, k.col(i) = res.msv.k;
        } y = x.transpose(), l = k.transpose();
    } else if (dim == 1) {
        x = res.msv.r, k = res.msv.k;
    }

    // calculate the potential values
    if (dim == 2) U = Expression(system.potential.at(0).at(0)).eval(x, y);
    if (dim == 1) U = Expression(system.potential.at(0).at(0)).eval(x);

    // define the initial wavefunction and normalize it
    if (dim == 2) {psi = Expression(opta.guess).eval(x, y); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dr);}
    if (dim == 1) {psi = Expression(opta.guess).eval(x); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dr);}

    // initialize the potential matrix
    res.msv.U = Matrix<>((int)std::pow(system.ngrid, dim), opta.real && !opta.spectrum.potential.empty() ? 2 : 1);

    // fill the ground state potential column
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
        // change the potential for te real time dynamics
        if (dim == 2) U = Expression(opta.spectrum.potential).eval(x, y).array() - (opta.spectrum.zpesub ? res.msv.opten(0) : 0);
        if (dim == 1) U = Expression(opta.spectrum.potential).eval(x).array() - (opta.spectrum.zpesub ? res.msv.opten(0) : 0);

        // fill the excited state potential column
        for (int i = 0; i < res.msv.U.rows(); i++) {
            res.msv.U(i, 1) = U(i % system.ngrid, i / system.ngrid);
        }

        // define the real time operators
        K = (-0.5 * I * (k.array().pow(2) + l.array().pow(2)) * opta.step / system.mass()).exp();
        R = (-0.5 * I * U.array() * opta.step).exp();
    }

    // loop over all states
    for (int i = 0; i < opta.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n  ITER         Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opta.real ? "(RT)" : "(IT)");

        // assign the guess as the optimized state to psi
        if (res.msv.optstates.size()) {psi = res.msv.optstates.at(i);}

        // calculate the total energy of the guess function
        Matrix<std::complex<double>> Ek = 0.5 * EigenConj(psi).array() * Numpy::FFT((k.array().pow(2) + l.array().pow(2)) * Numpy::FFT(psi).array(), 1).array();
        Matrix<std::complex<double>> Ep = EigenConj(psi).array() * U.array() * psi.array(); double E = (Ek / system.mass() + Ep).sum().real() * dr;

        // vector of state energy evolution, wfn vector and acf
        std::vector<Matrix<std::complex<double>>> wfns = {psi};
        Vector<std::complex<double>> acf(opta.iters + 1);
        std::vector<double> energies = {E}; acf(0) = 1;

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

            // push back the energies and wfn
            energies.push_back(E); if (opta.savewfn) wfns.push_back(psi);

            // add the point to acf
            if (opta.real && !opta.spectrum.potential.empty()) {
                acf(j) = ((EigenConj(wfns.at(0)).array() * psi.array()).sum() * dr);
            }

            // calculate the errors
            double Eerr = std::abs(E - Eprev), Derr = (psi.array().abs2() - Dprev.array()).abs2().sum();

            // print the iteration
            if (print) std::printf("%8d %20.14f %.2e %.2e\n", j, E, Eerr, Derr);
        }

        // assign the energy and wavefunction
        if (!opta.real) res.msv.optstates.at(i) = psi, res.msv.opten(i) = E;

        // block to calculate spectrum
        if (opta.real && !opta.spectrum.potential.empty()) {
            // append the ACF to the vector of ACFs
            res.msv.acfs.push_back(acf.array());

            // multiply the acf with the windowing function
            acf = acf.array() * Expression(opta.spectrum.window).eval(res.msv.t).array();

            // zero pad the acf
            acf.conservativeResize(acf.size() + opta.spectrum.zeropad); acf.tail(opta.spectrum.zeropad).setZero();

            // append and perform the Fourier transform
            res.msv.spectra.push_back(res.msv.f.array() * Numpy::HFFT(acf).array());

            // normalize the spectrum
            if (opta.spectrum.normalize) {
                res.msv.spectra.back() = res.msv.spectra.back().array() / res.msv.spectra.back().array().abs().maxCoeff();
            }
        }

        // save the state wavefunction
        if (opta.savewfn && dim == 2) ModelSystem::SaveWavefunction(ip + "/state" + std::to_string(i) + ".dat", x.col(0), y.row(0), wfns, energies);
        if (opta.savewfn && dim == 1) ModelSystem::SaveWavefunction(ip + "/state" + std::to_string(i) + ".dat", x, wfns, energies);
    }

    // print the newline
    if (print) std::printf("\n");

    // return the results
    return res;
}

Result ModelSolver::runnad(const ModelSystem& system, Result res, bool print) {
    // define the containers for energies, wfns, transforms, density matrices and coordinates and initialize wfn
    std::vector<double> energies; Matrix<std::complex<double>> psi(system.ngrid, 2); std::vector<Matrix<std::complex<double>>> psis;
    std::vector<Matrix<>> UT(system.ngrid, Matrix<>::Identity(system.potential.size(), system.potential.size()));
    res.msv.rho = Matrix<>(optn.iters, system.potential.size() * system.potential.size());
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid);

    // fill the real space and calculate dr
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } double dr = res.msv.r(1) - res.msv.r(0);

    // fill the momentum space
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / dr); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());

    // fill in the initial wavefunction and normalize it
    for (size_t i = 0; i < optn.guess.size(); i++) {
        psi.block(0, i, system.ngrid, 1) = (Expression(optn.guess.at(i)).eval(res.msv.r).array() * (I * optn.momentum * res.msv.r.array()).exp());
    } psi.block(0, 0, system.ngrid, 1).array() /= std::sqrt(psi.block(0, 0, system.ngrid, 1).array().abs2().sum() * dr);

    // create the potential
    std::vector<std::vector<Vector<std::complex<double>>>> V = {
        {Expression(system.potential.at(0).at(0)).eval(res.msv.r), Expression(system.potential.at(0).at(1)).eval(res.msv.r)},
        {Expression(system.potential.at(1).at(0)).eval(res.msv.r), Expression(system.potential.at(1).at(1)).eval(res.msv.r)}
    }; res.msv.U = Matrix<>(system.ngrid, system.potential.size());

    // subtract the complex absorbing potential
    V.at(0).at(0).array() -= (I * Expression(optn.cap).eval(res.msv.r)).array(), V.at(1).at(1).array() -= (I * Expression(optn.cap).eval(res.msv.r)).array();

    // loop to calculate the adiabatic transformation matrices
    for (int i = 0; i < system.ngrid; i++) {
        // define the diabatic potential at the current point
        Matrix<double> UD(V.size(), V.at(0).size());

        // fill the diabatic potential at the current point
        for (int j = 0; j < system.potential.size(); j++) {
            for (int k = 0; k < system.potential.at(j).size(); k++) {
                UD(j, k) = V.at(j).at(k)(i).real();
            }
        }

        // assign the diabatic potential to results
        res.msv.U.row(i) = UD.diagonal().transpose();
        
        // transform to adiabatic basis only if requested
        if (optn.adiabatic) {
            // diagonalize the diabatic potential to get the adiabatic potential
            Eigen::SelfAdjointEigenSolver<Matrix<>> eigensolver(UD); res.msv.U.row(i) = eigensolver.eigenvalues().transpose();

            // define the overlap of eigenvectors
            Vector<> overlap(system.potential.size());

            // calculate the overlap of eigenvectors
            for (int j = 0; j < system.potential.size(); j++) {
                overlap(j) = UT.at(std::max(i - 1, 0)).col(j).transpose() * eigensolver.eigenvectors().col(j); overlap(j) /= std::abs(overlap(j));
            }

            // maximize the overlap with the previuos transformation
            UT.at(i) = eigensolver.eigenvectors(); UT.at(i) = UT.at(i) * overlap.asDiagonal();
        }
    }

    // define the real and imaginary space operators
    std::vector<std::vector<Vector<std::complex<double>>>> K(2, std::vector<Vector<std::complex<double>>>(2)); auto R = K;

    // fill the kinetic propagator
    K.at(0).at(0) = (-0.5 * I * res.msv.k.array().pow(2) * optn.step / system.mass()).exp();
    K.at(1).at(1) = (-0.5 * I * res.msv.k.array().pow(2) * optn.step / system.mass()).exp();

    // fill the potential propagator
    Vector<std::complex<double>> D = 4 * V.at(1).at(0).array().abs().pow(2) + (V.at(0).at(0) - V.at(1).at(1)).array().pow(2);
    Vector<std::complex<double>> a = (-0.25 * I * (V.at(0).at(0) + V.at(1).at(1)) * optn.step).array().exp();
    Vector<std::complex<double>> b = (0.25 * D.array().sqrt() * optn.step).cos();
    Vector<std::complex<double>> c = I * (0.25 * D.array().sqrt() * optn.step).sin() / D.array().sqrt();
    R.at(0).at(0) = a.array() * (b.array() + c.array() * (V.at(1).at(1) - V.at(0).at(0)).array());
    R.at(0).at(1) = -2 * a.array() * c.array() * V.at(0).at(1).array();
    R.at(1).at(0) = -2 * a.array() * c.array() * V.at(0).at(1).array();
    R.at(1).at(1) = a.array() * (b.array() + c.array() * (V.at(0).at(0) - V.at(1).at(1)).array());

    // calculate the initial energy
    Vector<std::complex<double>> Ek00 = 0.5 * EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * Numpy::FFT(res.msv.k.array().pow(2) * Numpy::FFT(psi.col(0)).array(), 1).array();
    Vector<std::complex<double>> Ek11 = 0.5 * EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * Numpy::FFT(res.msv.k.array().pow(2) * Numpy::FFT(psi.col(1)).array(), 1).array();
    Vector<std::complex<double>> Ep00 = EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * V.at(0).at(0).array() * psi.col(0).array();
    Vector<std::complex<double>> Ep01 = EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * V.at(0).at(1).array() * psi.col(1).array();
    Vector<std::complex<double>> Ep10 = EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * V.at(1).at(0).array() * psi.col(0).array();
    Vector<std::complex<double>> Ep11 = EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * V.at(1).at(1).array() * psi.col(1).array();
    double E = ((Ek00 + Ek11) / system.mass() + Ep00 + Ep01 + Ep10 + Ep11).sum().real() * dr;

    // append the first wfn and energy
    energies.push_back(E), psis.push_back(psi);

    // print the iteration header
    if (print) std::printf(" ITER        Eel [Eh]         |dE|     |dD|     S1     S2\n");

    // propagate the wavefunction
    for (int i = 0; i < optn.iters; i++) {
        // save the previous values and temporary wavefunction
        Matrix<std::complex<double>> psip = psi, psit = psi; double Ep = E;

        // apply the first potential part of the propagator
        psit.col(0) = R.at(0).at(0).array() * psi.col(0).array() + R.at(0).at(1).array() * psi.col(1).array(); 
        psit.col(1) = R.at(1).at(0).array() * psi.col(0).array() + R.at(1).at(1).array() * psi.col(1).array();
        psi = psit;

        // apply the kinetic part of the propagator
        psit.col(0) = Numpy::FFT(K.at(0).at(0).array() * Numpy::FFT(psi.col(0)).array(), 1);
        psit.col(1) = Numpy::FFT(K.at(1).at(1).array() * Numpy::FFT(psi.col(1)).array(), 1);
        psi = psit;

        // apply the second potential part of the propagator
        psit.col(0) = R.at(0).at(0).array() * psi.col(0).array() + R.at(0).at(1).array() * psi.col(1).array();
        psit.col(1) = R.at(1).at(0).array() * psi.col(0).array() + R.at(1).at(1).array() * psi.col(1).array();
        psi = psit;

        // transform the temporary wavefunction to adiabatic basis
        for (int j = 0; j < psit.rows(); j++) psit.row(j) = UT.at(j).transpose() * psit.row(j).transpose();

        // calculate the adiabatic density matrix for the current time step
        Matrix<std::complex<double>> rho = psit.transpose() * psit.conjugate() * dr;

        // calculate the total energy
        Ek00 = 0.5 * EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * Numpy::FFT(res.msv.k.array().pow(2) * Numpy::FFT(psi.col(0)).array(), 1).array();
        Ek11 = 0.5 * EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * Numpy::FFT(res.msv.k.array().pow(2) * Numpy::FFT(psi.col(1)).array(), 1).array();
        Ep00 = EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * V.at(0).at(0).array() * psi.col(0).array();
        Ep01 = EigenConj<Vector<std::complex<double>>>(psi.col(0)).array() * V.at(0).at(1).array() * psi.col(1).array();
        Ep10 = EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * V.at(1).at(0).array() * psi.col(0).array();
        Ep11 = EigenConj<Vector<std::complex<double>>>(psi.col(1)).array() * V.at(1).at(1).array() * psi.col(1).array();
        E = ((Ek00 + Ek11) / system.mass() + Ep00 + Ep01 + Ep10 + Ep11).sum().real() * dr;

        // calculate the errors and assign the previous values
        double Eerr = std::abs(E - Ep), Derr = (psi.array().abs2() - psip.array().abs2()).abs2().sum();

        // fill the density matrix container
        res.msv.rho(i, 1) = std::abs(rho(0, 0));
        res.msv.rho(i, 2) = std::abs(rho(1, 1));
        res.msv.rho(i, 3) = std::abs(rho(0, 1));
        res.msv.rho(i, 0) = i * optn.step;

        // append the temporary wfn in adiabatic basis and energy
        energies.push_back(E); if (optn.savewfn) psis.push_back(psit);

        // print the iteration
        if (print) std::printf("%6d %20.14f %.2e %.2e %6.4f %6.4f\n", i + 1, E, Eerr, Derr, res.msv.rho(i, 1), res.msv.rho(i, 2));
    }

    // save the wfn dynamics and potential in adiabatic basis
    for (size_t i = 0; i < system.potential.size() && optn.savewfn; i++) {
        // extract the wavefunction in diabatic basis
        std::vector<Matrix<std::complex<double>>> state; for (auto psi : psis) state.push_back(psi.col(i));

        // save the wavefunction
        ModelSystem::SaveWavefunction(ip + "/state" + std::to_string(i) + ".dat", res.msv.r, state, energies);
    }

    // print the newline
    if (print) std::printf("\n");

    // return
    return res;
}

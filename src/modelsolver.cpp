#include "modelsolver.h"
#include "modelsystem.h"

#define I std::complex<double>(0, 1)

Result ModelSolver::run(const ModelSystem& system, Result res, bool print) {
    // define the space, momentum, time and frequency grids
    res.msv.r = Vector<>(system.ngrid), res.msv.k = Vector<>(system.ngrid), res.msv.t = Vector<>(opta.iters + 1), res.msv.f = Vector<>(opta.iters + opta.spectrum.zeropad + 1);

    // fill the real and time domain
    for (int i = 0; i < system.ngrid; i++) {
        res.msv.r(i) = system.limits.at(0) + (system.limits.at(1) - system.limits.at(0)) * i / (system.ngrid - 1);
    } for (int i = 0; i <= opta.iters; i++) res.msv.t(i) = i * opta.step;

    // fill the momentum and frequency domain
    res.msv.k.fill(2 * M_PI / res.msv.k.size() / (res.msv.r(1) - res.msv.r(0))); for (int i = 0; i < res.msv.k.size(); i++) res.msv.k(i) *= i - (i < res.msv.k.size() / 2 ? 0 : res.msv.k.size());
    res.msv.f.fill(2 * M_PI / res.msv.f.size() / (res.msv.t(1) - res.msv.t(0))); for (int i = 0; i < res.msv.f.size(); i++) res.msv.f(i) *= i - (i < res.msv.f.size() / 2 ? 0 : res.msv.f.size());

    // calculate the space element
    res.msv.dr = std::pow(res.msv.r(1) - res.msv.r(0), system.vars().size());

    // initialize the potential matrix
    res.msv.U = Matrix<>((int)std::pow(system.ngrid, system.vars().size()), system.potential.size() + (opta.real && !opta.spectrum.potential.empty() ? 1 : 0));

    // run the calculation
    if (system.potential.size() == 1) return runad(system, res, print);
    else return runnad(system, res, print);
}

Result ModelSolver::runad(const ModelSystem& system, Result res, bool print) {
    // run the classical dynamics if requested
    if (!optd.gradient.empty()) return runcd(system, res, print);

    // define the wavefunction, potential and independent variables container
    Matrix<std::complex<double>> psi; Matrix<> x, y, k, l; RealMatrixOfMatrices V(1, 1);

    // initialize the 2-dimensional variables
    if (system.vars().size() == 2) {
        x = Matrix<>::Zero(system.ngrid, system.ngrid), y = Matrix<>::Zero(system.ngrid, system.ngrid);
        k = Matrix<>::Zero(system.ngrid, system.ngrid), l = Matrix<>::Zero(system.ngrid, system.ngrid);
    } else if (system.vars().size() == 1) {
        x = Vector<>::Zero(system.ngrid), y = Vector<>::Zero(system.ngrid);
        k = Vector<>::Zero(system.ngrid), l = Vector<>::Zero(system.ngrid);
    }

    // fill the independent variables
    if (system.vars().size() == 2) {
        for (int i = 0; i < system.ngrid; i++) {
            x.col(i) = res.msv.r, k.col(i) = res.msv.k;
        } y = x.transpose(), l = k.transpose();
    } else if (system.vars().size() == 1) {
        x = res.msv.r, k = res.msv.k;
    }

    // calculate the potential values
    if (system.vars().size() == 2) V(0, 0) = Expression(system.potential.at(0).at(0), system.vars()).eval(x, y);
    if (system.vars().size() == 1) V(0, 0) = Expression(system.potential.at(0).at(0), system.vars()).eval(x);

    // define the initial wavefunction and normalize it
    if (system.vars().size() == 2) {psi = Expression(opta.guess, system.vars()).eval(x, y); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * res.msv.dr);}
    if (system.vars().size() == 1) {psi = Expression(opta.guess, system.vars()).eval(x); psi = psi.array() / std::sqrt(psi.array().abs2().sum() * res.msv.dr);}

    // fill the ground state potential column
    for (int i = 0; i < res.msv.U.rows(); i++) {
        res.msv.U(i, 0) = V(0, 0)(i % system.ngrid, i / system.ngrid);
    }

    // define the vector of optimal wavefunction and energies if imaginary time is selected
    if (!opta.real) res.msv.optstates = std::vector<Matrix<std::complex<double>>>(opta.nstate, psi), res.msv.opten = Vector<>(opta.nstate); 

    // define the imaginary time operators
    auto[R, K] = Propagator<1>::Get(system, V.complex(), k.array().pow(2) + l.array().pow(2), opta.step, false);

    // real time dynamics operators
    if (opta.real) {
        // change the potential for the real time dynamics
        if (system.vars().size() == 2) V(0, 0) = Expression(opta.spectrum.potential, system.vars()).eval(x, y).array() - (opta.spectrum.zpesub ? res.msv.opten(0) : 0);
        if (system.vars().size() == 1) V(0, 0) = Expression(opta.spectrum.potential, system.vars()).eval(x).array() - (opta.spectrum.zpesub ? res.msv.opten(0) : 0);

        // fill the excited state potential column
        for (int i = 0; i < res.msv.U.rows(); i++) {
            res.msv.U(i, 1) = V(0, 0)(i % system.ngrid, i / system.ngrid);
        }

        // define the real time operators
        std::tie(R, K) = Propagator<1>::Get(system, V.complex(), k.array().pow(2) + l.array().pow(2), opta.step, true);
    }

    // loop over all states
    for (int i = 0; i < opta.nstate; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i %s\n  ITER         Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i, opta.real ? "(RT)" : "(IT)");

        // assign the guess as the optimized state to psi
        if (res.msv.optstates.size()) {psi = res.msv.optstates.at(i);}

        // calculate the total energy of the guess function
        double E = Propagator<1>::Energy(system, V, res.msv.r, k.array().pow(2) + l.array().pow(2), psi);

        // vector of state energy evolution, wfn vector and acf
        std::vector<Matrix<std::complex<double>>> wfns = {psi};
        Vector<std::complex<double>> acf(opta.iters + 1);
        std::vector<double> energies = {E}; acf(0) = 1;

        // propagate the current state
        for (int j = 1; j <= opta.iters; j++) {
            // save the previous values
            Matrix<std::complex<double>> Dprev = psi.array().abs2(); double Eprev = E;

            // propagate the wavefunction
            psi = Propagator<1>::Propagate(R, K, psi);

            // subtract lower eigenstates
            for (int k = 0; k < i && !opta.real; k++) {
                psi = psi.array() - (EigenConj(res.msv.optstates.at(k)).array() * psi.array()).sum() * res.msv.dr * res.msv.optstates.at(k).array();
            }

            // normalize the wavefunction
            psi = psi.array() / std::sqrt(psi.array().abs2().sum() * res.msv.dr);

            // calculate the total energy
            E = Propagator<1>::Energy(system, V, res.msv.r, k.array().pow(2) + l.array().pow(2), psi);

            // push back the energies and wfn
            energies.push_back(E); if (opta.savewfn) wfns.push_back(psi);

            // add the point to acf
            if (opta.real && !opta.spectrum.potential.empty()) {
                acf(j) = ((EigenConj(wfns.at(0)).array() * psi.array()).sum() * res.msv.dr);
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
            acf = acf.array() * Expression(opta.spectrum.window, system.vars()).eval(res.msv.t).array();

            // zero pad the acf
            acf.conservativeResize(acf.size() + opta.spectrum.zeropad); acf.tail(opta.spectrum.zeropad).setZero();

            // append and perform the Fourier transform
            res.msv.spectra.push_back(Numpy::HFFT(acf).array());

            // normalize the spectrum
            if (opta.spectrum.normalize) {
                res.msv.spectra.back() = res.msv.spectra.back().array() / res.msv.spectra.back().array().abs().maxCoeff();
            }
        }

        // save the state wavefunction
        if (opta.savewfn && system.vars().size() == 2) ModelSystem::SaveWavefunction(ip / ("state" + std::to_string(i + 1) + ".dat"), x.col(0), y.row(0), wfns, energies);
        if (opta.savewfn && system.vars().size() == 1) ModelSystem::SaveWavefunction(ip / ("state" + std::to_string(i + 1) + ".dat"), x, wfns, energies);
    }

    // print the newline
    if (print) std::printf("\n");

    // return the results
    return res;
}

Result ModelSolver::runnad(const ModelSystem& system, Result res, bool print) {
    // run the classical dynamics if requested
    if (!optd.gradient.empty()) return runcd(system, res, print);

    // define the containers for energies, wfns, transforms, density matrices and coordinates and initialize wfn
    std::vector<double> energies; Matrix<std::complex<double>> psi(system.ngrid, 2); std::vector<Matrix<std::complex<double>>> psis;
    std::vector<Matrix<>> UT(system.ngrid, Matrix<>::Identity(system.potential.size(), system.potential.size()));
    res.msv.rho = Matrix<>(optn.iters, system.potential.size() * system.potential.size());

    // fill in the initial wavefunction and normalize it
    for (size_t i = 0; i < optn.guess.size(); i++) {
        psi.block(0, i, system.ngrid, 1) = (Expression(optn.guess.at(i), system.vars()).eval(res.msv.r).array() * (I * optn.momentum * res.msv.r.array()).exp());
    } psi.block(0, 0, system.ngrid, 1).array() /= std::sqrt(psi.block(0, 0, system.ngrid, 1).array().abs2().sum() * res.msv.dr);

    // create and fill the potential
    ComplexMatrixOfMatrices V(2, 2);
    V(0, 0) = Expression(system.potential.at(0).at(0), system.vars()).eval(res.msv.r);
    V(0, 1) = Expression(system.potential.at(0).at(1), system.vars()).eval(res.msv.r);
    V(1, 0) = Expression(system.potential.at(1).at(0), system.vars()).eval(res.msv.r);
    V(1, 1) = Expression(system.potential.at(1).at(1), system.vars()).eval(res.msv.r);

    // subtract the complex absorbing potential
    V(0, 0).array() -= (I * Expression(optn.cap, system.vars()).eval(res.msv.r)).array(), V(1, 1).array() -= (I * Expression(optn.cap, system.vars()).eval(res.msv.r)).array();

    // loop to calculate the adiabatic transformation matrices
    for (int i = 0; i < system.ngrid; i++) {
        // define the diabatic potential at the current point
        Matrix<double> UD(system.potential.size(), system.potential.size());

        // fill the diabatic potential at the current point
        for (size_t j = 0; j < system.potential.size(); j++) {
            for (size_t k = 0; k < system.potential.at(j).size(); k++) {
                UD(j, k) = V(j, k)(i).real();
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
            for (size_t j = 0; j < system.potential.size(); j++) {
                overlap(j) = UT.at(std::max(i - 1, 0)).col(j).transpose() * eigensolver.eigenvectors().col(j); overlap(j) /= std::abs(overlap(j));
            }

            // maximize the overlap with the previuos transformation
            UT.at(i) = eigensolver.eigenvectors(); UT.at(i) = UT.at(i) * overlap.asDiagonal();
        }
    }

    // define the real and imaginary space operators
    auto[R, K] = Propagator<2>::Get(system, V.complex(), res.msv.k.array().pow(2), optn.step, true);

    // calculate the initial energy
    double E = Propagator<2>::Energy(system, V.real(), res.msv.r, res.msv.k.array().pow(2), psi);

    // append the first wfn and energy
    energies.push_back(E), psis.push_back(psi);

    // print the iteration header
    if (print) std::printf(" ITER        Eel [Eh]         |dE|     |dD|     S1     S2\n");

    // propagate the wavefunction
    for (int i = 0; i < optn.iters; i++) {
        // save the previous values and temporary wavefunction
        Matrix<std::complex<double>> psip = psi; double Ep = E;

        // apply the propagator
        psi = Propagator<2>::Propagate(R, K, psi);

        // copy the wavefunction and transform it to the adiabatic basis if the matrices UT were calculated
        Matrix<std::complex<double>> psitemp = psi; for (int j = 0; j < psitemp.rows(); j++) psitemp.row(j) = UT.at(j).transpose() * psitemp.row(j).transpose();

        // calculate the adiabatic density matrix for the current time step
        Matrix<std::complex<double>> rho = psitemp.transpose() * psitemp.conjugate() * res.msv.dr;

        // calculate the total energy
        E = Propagator<2>::Energy(system, V.real(), res.msv.r, res.msv.k.array().pow(2), psi);

        // calculate the errors and assign the previous values
        double Eerr = std::abs(E - Ep), Derr = (psi.array().abs2() - psip.array().abs2()).abs2().sum();

        // fill the density matrix container
        res.msv.rho(i, 1) = std::abs(rho(0, 0));
        res.msv.rho(i, 2) = std::abs(rho(1, 1));
        res.msv.rho(i, 3) = std::abs(rho(0, 1));
        res.msv.rho(i, 0) = i * optn.step;

        // append the temporary wfn in adiabatic basis and energy
        energies.push_back(E); if (optn.savewfn) psis.push_back(psitemp);

        // print the iteration
        if (print) std::printf("%6d %20.14f %.2e %.2e %6.4f %6.4f\n", i + 1, E, Eerr, Derr, res.msv.rho(i, 1), res.msv.rho(i, 2));
    }

    // save the wfn dynamics and potential in chosen basis
    for (size_t i = 0; i < system.potential.size() && optn.savewfn; i++) {
        // extract the wavefunction in chosen basis
        std::vector<Matrix<std::complex<double>>> state; for (auto psi : psis) state.push_back(psi.col(i));

        // save the wavefunction
        ModelSystem::SaveWavefunction(ip / ("state" + std::to_string(i + 1) + ".dat"), res.msv.r, state, energies);
    }

    // print the newline
    if (print) std::printf("\n");

    // return
    return res;
}

Result ModelSolver::runcd(const ModelSystem& system, Result res, bool print) {
    // print the header
    if (print) std::printf(" TRAJ   ITER  TIME [fs] STATE PROB    POT [Eh]       KIN [Eh]        E [Eh]      |GRAD|      TIME\n");

    // define energy, gradient and coupling functions and define the array of state populations
    std::vector<Expression> energy, coupling, grad; std::vector<Vector<int>> states(optd.trajs);

    // fill the gradient and energy expressions
    for (size_t i = 0; i < system.potential.size(); i++) {
        energy.push_back(Expression(system.potential.at(i).at(i), system.vars()));
        grad.push_back(Expression(optd.gradient.at(i), system.vars()));
        res.msv.U.col(i) = energy.at(i).eval(res.msv.r);
    }

    // fill the coupling expressions
    for (size_t i = 0; i < system.potential.size(); i++) {
        for (size_t j = 0; j < system.potential.size(); j++) {
            if (i != j) coupling.push_back(Expression(system.potential.at(i).at(j), system.vars()));
        }
    }

    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < optd.trajs; i++) {
        // random number generators
        std::mt19937 mt(optd.seed + i); std::uniform_real_distribution<double> dist(0, 1);
    
        // create the containers for the position, velocity, acceleration and mass
        Matrix<> r(optd.iters + 1, system.vars().size()); Vector<> v(system.vars().size()), a(system.vars().size()), m(system.vars().size()); a.fill(0), m.fill(system.mass());

        // fill the initial position and momentum
        std::normal_distribution<double> positiondist(optd.position.at(0), 0.5);
        std::normal_distribution<double> momentumdist(optd.momentum.at(0), 1.0);
        r.row(0)(0) = positiondist(mt), v(0) = momentumdist(mt) / system.mass();

        // define the initial state, state container and start the timer
        Vector<int> state(optd.iters + 1); state(0) = optd.state - 1; auto start = Timer::Now();

        // calculate the force and potential with kinetic energy
        Vector<> F = -grad.at(state(0)).eval(r.row(0)); double Epot = energy.at(state(0)).get(r.row(0)), Ekin = 0.5 * (m.array() * v.array() * v.array()).sum();

        // define the energy difference container and fill the zeroth iteration
        std::vector<double> Ediff; if (system.potential.size() > 1) Ediff.push_back(energy.at(1).get(r.row(0)) - energy.at(0).get(r.row(0)));

        // print the zeroth iteration
        if (print) std::printf("%6d %6d %9.4f %5d %1.2f %14.8f %14.8f %14.8f %.2e %s\n", i + 1, 0, 0.0, state(0) + 1, 0.0, Epot, Ekin, Epot + Ekin, F.norm(), Timer::Format(Timer::Elapsed(start)).c_str());

        for (int j = 0; j < optd.iters; j++) {
            // start the timer and store the previous v and a
            start = Timer::Now(); Vector<> vp = v, ap = a;

            // calculate the velocity and accceleration
            a = F.array() / m.array(); v = vp + 0.5 * (ap + a) * optd.step;

            // move the system
            r.row(j + 1) = r.row(j) + optd.step * (v + 0.5 * a * optd.step);

            // fill the state
            state(j + 1) = state(j);

            // calculate the potential and kinetic energy
            double Epot = energy.at(state(j + 1)).get(r.row(j + 1)), Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(), P = 0;

            if (system.potential.size() > 1) {
                // append the energy and calculate the derivative of the energy difference
                Ediff.push_back(energy.at(1).get(r.row(j + 1)) - energy.at(0).get(r.row(j + 1))); double dEdiff = (Ediff.at(j + 1) - Ediff.at(j)) / optd.step;

                // calculate the probability of state change according to the Landau-Zener formula
                double gamma = std::pow(coupling.at(0).get(r.row(j + 1)), 2) / std::abs(dEdiff); P = 1 - std::exp(-2 * M_PI * gamma);

                // change the state if the jump is accepted
                if (Ediff.at(j) * Ediff.at(j + 1) < 0 && dist(mt) < P) {
                    state(j + 1) = state(j + 1) == 1 ? 0 : 1; v(0) = std::sqrt(v(0)*v(0) - (state(j + 1) - state(j)) * 2 * Ediff.at(j + 1) / system.mass());
                }
            }

            // calculate the force
            F = -grad.at(state(j + 1)).eval(r.row(j + 1));

            // print the iteration
            if (print) std::printf("%6d %6d %9.4f %5d %1.2f %14.8f %14.8f %14.8f %.2e %s\n", i + 1, j + 1, AU2FS * optd.step * (j + 1), state(j) + 1, P, Epot, Ekin, Epot + Ekin, F.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
        }

        // save the states and define the density diagonal
        states.at(i) = state;

        // save the matrix of states and coordinates
        Matrix<> SR(r.rows(), r.cols() + 1); SR << (state.array() + 1).cast<double>(), r; EigenWrite(ip / ("trajectory" + std::to_string(i + 1) + ".mat"), SR);
    }

    if (system.potential.size() > 1) {
        // define the density matrix diagonal array
        Matrix<> rho = Matrix<>(optd.iters + 1, system.potential.size() + 1);

        // fill the density matrix diagonal
        for (int j = 0; j < optd.iters + 1; j++) {
            for (const auto& state : states) {
                rho(j, 1) += state(j) == 0 ? 1.0 / states.size() : 0;
            } rho(j, 0) = j * optd.step, rho(j, 2) = 1 - rho(j, 1);
        } EigenWrite(ip / ("rho.mat"), rho);
    }


    // print the newline
    if (print) std::printf("\n");

    // return the results
    return res;
}

#include "qdyn.h"

void Acorn::QDYN::run(const Options& opt, std::vector<timepoint>& timers) {
    // start the timer for the potential and wavefunction loading
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of potential and wavefunction loading
    std::cout << "READING POTENTIAL AND INITIAL WAVEFUNCTION: " << std::flush;

    // load the potential and initial wavefunction
    Eigen::MatrixXd Ud = torch::ReadMatrix("U_DIA.mat"); std::vector<Acorn::QDYN::Wavefunction> states(opt.optstates, Acorn::QDYN::Wavefunction(torch::ReadMatrix("PSI_DIA_GUESS.mat").rightCols(2 * (int)std::sqrt(Ud.cols() - opt.dim)), Ud.leftCols(opt.dim), opt.mass, opt.momentum).normalized());

    // print the time for potential and wavefunction loading
    std::cout << eltime(timers.at(1)) << std::endl;

    // extract the number of states
    int nstate = states.at(0).get().cols(), nstatessq = nstate * nstate;

    // remove the first column from matrix (the independent variable column)
    Ud.block(0, 0, Ud.rows(), Ud.cols() - opt.dim) = Ud.rightCols(Ud.cols() - opt.dim); Ud.conservativeResize(Ud.rows(), Ud.cols() - opt.dim);

    // define and initialize the adiabatic variables
    Eigen::MatrixXd Ua; std::vector<Eigen::MatrixXd> UT; if (opt.adiabatic) {
        Ua = Eigen::MatrixXd(Ud.rows(), nstate + opt.dim); UT = std::vector<Eigen::MatrixXd>(Ud.rows(), Eigen::MatrixXd::Identity(nstate, nstate));
    }

    // diagonalize the potential at each point
    for (int i = 0; i < Ud.rows() && opt.adiabatic; i++) {

        // initialize the diabatic potential at the current point
        Eigen::MatrixXd UD = Ud.row(i).reshaped(nstate, nstate).transpose();

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(UD); Eigen::MatrixXd C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(nstate, 1);

        // fill the diabatic matrix
        Ua(i, 0) = states.at(0).getr()(i); Ua.row(i).rightCols(nstate) = eps.transpose();

        // calculate the overlap of eigenvectors
        for (int j = 0; j < nstate; j++) {
            overlap(j) = UT.at(std::max(i - 1, 0)).col(j).transpose() * C.col(j); overlap(j) = overlap(j) < 0 ? -1 : 1;
        }

        // maximize the overlap with the previuos transformation
        UT.at(i) = C * overlap.asDiagonal();
    }

    // save the adiabatic potential
    if (opt.adiabatic) {
    
        // start the timer for saving the adiabatic potential
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of saving the adiabatic potential
        std::cout << "\nSAVING THE ADIABATIC POTENTIAL: " << std::flush;

        // save the adiabatic potential
        torch::WriteMatrix("U_ADIA.mat", Ua);

        // print the time for saving the adiabatic potential
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // define the real time kinetic and potential operators in diabatic basis
    auto[R, K] = states.at(0).propagator(Ud, std::complex<double>(opt.imaginary, !opt.imaginary), opt.step);

    // define the energy, wfn and acf vector
    Eigen::VectorXd eps(opt.optstates); Eigen::MatrixXd acf(opt.imaginary ? 0 : opt.iters + 1, opt.imaginary ? 0 : 3);

    // for every orthogonal state
    for (int i = 0; i < opt.optstates; i++) {

        // set the initial wfn, energy and acf value
        Acorn::QDYN::Wavefunction wfnd = states.at(i); double E = wfnd.energy(Ud); if (!opt.imaginary) acf(0, 0) = 0, acf(0, 1) = 1, acf(0, 2) = 0;

        // define the diabatic and adiabatic wavefunction and density matrix containers
        Acorn::QDYN::Wavefunction wfna; Eigen::MatrixXd wfndt, wfnat, Pd, Pa, rho; if (nstate > 1) Pd = Eigen::MatrixXd(opt.iters + 1, nstatessq + 1), Pa = Eigen::MatrixXd(opt.iters + 1, nstatessq + 1), rho = Eigen::MatrixXd(nstate, nstate);

        // transform the wfn to the adiabatic basis
        if (opt.adiabatic) wfna = wfnd.adiabatize(UT);

        // save the diabatic wavefunction and fill the independent variable columns
        if (opt.savewfn) {
            wfndt = Eigen::MatrixXd(Ud.rows(), 2 * nstate * (opt.iters + 1) + opt.dim);
            for (int i = 0; i < nstate; i++) {
                wfndt.col(2 * i + opt.dim) = opt.factor * wfnd.get().col(i).real(), wfndt.col(2 * i + opt.dim + 1) = opt.factor * wfnd.get().col(i).imag();
                if (opt.align) {
                    wfndt.col(2 * i + opt.dim + 0) = wfndt.col(2 * i + opt.dim + 0).array() + Ud.col(i * (nstate + 1)).array();
                    wfndt.col(2 * i + opt.dim + 1) = wfndt.col(2 * i + opt.dim + 1).array() + Ud.col(i * (nstate + 1)).array();
                }
            } wfndt.leftCols(opt.dim) = wfnd.getr();
        }

        // save the adiabatic wavefunction and fill the independent variable columns
        if (opt.adiabatic && opt.savewfn) {
            wfnat = Eigen::MatrixXd(Ua.rows(), 2 * nstate * (opt.iters + 1) + opt.dim);
            for (int i = 0; i < nstate; i++) {
                wfnat.col(2 * i + opt.dim) = opt.factor * wfna.get().col(i).real(), wfnat.col(2 * i + opt.dim + 1) = opt.factor * wfna.get().col(i).imag();
                if (opt.align) {
                    wfnat.col(2 * i + opt.dim + 0) = wfnat.col(2 * i + opt.dim + 0).array() + Ua.col(1 + i).array();
                    wfnat.col(2 * i + opt.dim + 1) = wfnat.col(2 * i + opt.dim + 1).array() + Ua.col(1 + i).array();
                }
            } wfnat.leftCols(opt.dim) = wfna.getr();
        }

        // calculate the diabatic and adiabatic density matrices and save them
        if (nstate > 1) {
            if (opt.adiabatic) {
                rho = wfna.density(); Pa.row(0)(0) = 0, Pa.row(0).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);
            } rho = wfnd.density(); Pd.row(0)(0) = 0, Pd.row(0).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);
        }

        // print the header
        std::printf("\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIME");

        // perform the propagation
        for (int j = 0; j < opt.iters; j++) {

            // reset the timer, save the previous energy and propagate
            timers.at(1) = std::chrono::high_resolution_clock().now(); double Ep = E; wfnd = wfnd.propagate(R, K);

            // subrtract the lower eigenstates in ITP
            for (int k = 0; k < i && opt.imaginary; k++) wfnd = wfnd - states.at(k) * states.at(k).overlap(wfnd);

            // normalize the wavefunction if ITP and calculate the energy
            if (opt.imaginary) {wfnd = wfnd.normalized();} E = wfnd.energy(Ud);

            // assign the ACF value
            std::complex<double> overlap; if (!opt.imaginary) {overlap = states.at(i).overlap(wfnd); acf(j + 1, 0) = (j + 1) * opt.step; acf(j + 1, 1) = overlap.real(); acf(j + 1, 2) = overlap.imag();}

            // transform the wfn to the adiabatic basis and save the adiabatic density matrix
            if (opt.adiabatic) {
                wfna = wfnd.adiabatize(UT); if (nstate > 1) {rho = wfna.density(); Pa(j + 1, 0) = (j + 1) * opt.step, Pa.row(j + 1).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);}
            }

            // save the diabatic density matrix
            if (nstate > 1) {rho = wfnd.density(); Pd(j + 1, 0) = (j + 1) * opt.step, Pd.row(j + 1).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);}

            // save the diabatic wavefunction to the wavefunction containers
            if (opt.savewfn) {
                for (int k = 0; k < nstate; k++) {
                    wfndt.col(2 * nstate * (j + 1) + 2 * k + opt.dim) = opt.factor * wfnd.get().col(k).real(), wfndt.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 1) = opt.factor * wfnd.get().col(k).imag();
                    if (opt.align) {
                        wfndt.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 0) = wfndt.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 0).array() + Ud.col(k * (nstate + 1)).array();
                        wfndt.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 1) = wfndt.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 1).array() + Ud.col(k * (nstate + 1)).array();
                    }
                }
            }

            // save the adiabatic wavefunction to the wavefunction containers
            if (opt.adiabatic && opt.savewfn) {
                for (int k = 0; k < nstate; k++) {
                    wfnat.col(2 * nstate * (j + 1) + 2 * k + opt.dim) = opt.factor * wfna.get().col(k).real(), wfnat.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 1) = opt.factor * wfna.get().col(k).imag();
                    if (opt.align) {
                        wfnat.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 0) = wfnat.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 0).array() + Ua.col(1 + k).array();
                        wfnat.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 1) = wfnat.col(2 * nstate * (j + 1) + 2 * k + opt.dim + 1).array() + Ua.col(1 + k).array();
                    }
                }
            }

            // print the iteration info
            std::printf("%6d %20.14f %.2e %s\n", j + 1, E, std::abs(E - Ep), eltime(timers.at(1)).c_str());
        }

        // assign results and print a new line
        states.at(i) = wfnd; eps(i) = E; std::cout << std::endl;

        // define the spectrum matrix
        Eigen::MatrixXd spectrum(opt.imaginary ? 0 : opt.iters + 1, opt.imaginary ? 0 : 2);

        // calculate the spectrum
        if (!opt.imaginary) {

            // perform the fourier transform of the ACF
            spectrum.col(1) = Acorn::FourierTransform::C2RFFT(acf.col(1) + std::complex<double>(0, 1) * acf.col(2)).array().abs();

            // fill the frequency column of the spectrum
            spectrum.col(0).fill(2 * M_PI / (opt.iters + 1) / opt.step); for (int i = 0; i < opt.iters + 1; i++) spectrum.col(0)(i) *= i - (i < (opt.iters + 1) / 2 ? 0 : (opt.iters + 1));
        }

        // print the populations
        for (int i = 0; i < nstate && nstate > 1; i++) {
            std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(nstatessq).reshaped(nstate, nstate).transpose()(i, i));
        } if (nstate > 1) std::cout << std::endl;

        // print the populations
        for (int i = 0; i < nstate && nstate > 1 && opt.adiabatic; i++) {
            std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(nstatessq).reshaped(nstate, nstate).transpose()(i, i));
        } if (nstate > 1 && opt.adiabatic) std::cout << std::endl;

        // start the timer for writing the results
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of the writing of final results
        std::cout << "WAVEFUNCTIONS, DENSITY MATRICES, ACF AND SPECTRUM WRITING: " << std::flush;

        // save the resulting data
        if (                 opt.savewfn) torch::WriteMatrix("PSI_DIA_"      + std::to_string(i) + ".mat", wfndt   );
        if (nstate > 1                  ) torch::WriteMatrix("P_DIA_"        + std::to_string(i) + ".mat", Pd      );
        if (opt.adiabatic && opt.savewfn) torch::WriteMatrix("PSI_ADIA_"     + std::to_string(i) + ".mat", wfnat   );
        if (opt.adiabatic               ) torch::WriteMatrix("P_ADIA_"       + std::to_string(i) + ".mat", Pa      );
        if (!opt.imaginary              ) torch::WriteMatrix("ACF_DIA_"      + std::to_string(i) + ".mat", acf     );
        if (!opt.imaginary              ) torch::WriteMatrix("SPECTRUM_DIA_" + std::to_string(i) + ".mat", spectrum);

        // print the time for writing the results
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // print the state energies
    std::cout << std::endl; for (int i = 0; i < opt.optstates; i++) std::printf("FINAL WFN %02d ENERGY: %.14f\n", i, eps(i));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << eltime(timers.at(0)) << std::endl;
}

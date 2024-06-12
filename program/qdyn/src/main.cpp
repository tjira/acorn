#include "fourier.h"
#include "timer.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Nonadiabatic Quantum Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-o", "--optimize").help("-- Perform the dynamics in imaginary time for several orthogonal states.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("--adiabatic").help("-- Enable transform to adiabatic basis.").default_value(false).implicit_value(true);
    program.add_argument("--savewfn").help("-- Save the time evolution of the wavefunction.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters
    int iters = program.get<int>("-i"); double mass = program.get<double>("-m"), step = program.get<double>("-s"), momentum = program.get<double>("-p");

    // extract boolean flags
    bool adiabatic = program.get<bool>("--adiabatic"), imaginary = program.is_used("--optimize"), savewfn = program.get<bool>("--savewfn");

    // load the potential in the diabatic basis
    MEASURE("INITIAL WAVEFUNCTION AND POTENTIAL IN DIABATIC BASIS READING: ",
        Matrix Ud = Eigen::LoadMatrix("U_DIA.mat"); Wavefunction wfnd(Eigen::LoadMatrix("PSI_DIA_GUESS.mat").rightCols(2 * (int)std::sqrt(Ud.cols() - 1)), Ud.leftCols(1), mass, momentum);
    )

    // extract the number of states
    int nstate = wfnd.states(), nstatessq = nstate * nstate;

    // normalize the wfn, remove first column from matrix (the independent variable column)
    wfnd = wfnd.normalized(); Ud.block(0, 0, Ud.rows(), Ud.cols() - 1) = Ud.rightCols(Ud.cols() - 1); Ud.conservativeResize(Ud.rows(), Ud.cols() - 1);

    // define and initialize the adiabatic variables
    Matrix Ua; std::vector<Matrix> UT; if (adiabatic) {
        Ua = Matrix(Ud.rows(), nstate + 1); UT = std::vector<Matrix>(Ud.rows(), Matrix::Identity(nstate, nstate));
    }

    // diagonalize the potential at each point
    for (int i = 0; i < Ud.rows() && adiabatic; i++) {

        // initialize the diabatic potential at the current point
        Matrix UD = Ud.row(i).reshaped(nstate, nstate).transpose();

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Matrix> solver(UD); Matrix C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(nstate, 1);

        // fill the diabatic matrix
        Ua(i, 0) = wfnd.getr()(i); Ua.row(i).rightCols(nstate) = eps.transpose();

        // calculate the overlap of eigenvectors
        for (int j = 0; j < nstate; j++) {
            overlap(j) = UT.at(std::max(i - 1, 0)).col(j).transpose() * C.col(j); overlap(j) = overlap(j) < 0 ? -1 : 1;
        }

        // maximize the overlap with the previuos transformation
        UT.at(i) = C * overlap.asDiagonal();
    }

    // define the real time kinetic and potential operators in diabatic basis
    auto[R, K] = wfnd.propagator(Ud, std::complex<double>(imaginary, !imaginary), step);

    // define the energy, wfn and acf vector
    Vector eps(program.get<int>("-o")); Matrix acf(iters + 1, 3); std::vector<Wavefunction> states(program.get<int>("-o"), wfnd);

    // for every orthogonal state
    for (int i = 0; i < program.get<int>("-o"); i++) {

        // set the initial wfn, energy and acf value
        wfnd = states.at(i); double E = wfnd.energy(Ud); acf(0, 0) = 0, acf(0, 1) = 1, acf(0, 2) = 0;

        Wavefunction wfna; Matrix wfndt, wfnat, Pd(iters + 1, nstatessq + 1), Pa(iters + 1, nstatessq + 1), rho(nstate, nstate);

        if (adiabatic) wfna = wfnd.adiabatize(UT);

        if (savewfn) {
            wfndt = Matrix(Ud.rows(), 2 * nstate * (iters + 1) + 1);
            for (int i = 0; i < nstate; i++) {
                wfndt.col(2 * i + 1) = wfnd.get().col(i).real(), wfndt.col(2 * i + 2) = wfnd.get().col(i).imag();
            } wfndt.col(0) = wfnd.getr();
        }

        // adiabatize the guess wavefunction and fill the independent variable columns
        if (adiabatic && savewfn) {
            wfnat = Matrix(Ua.rows(), 2 * nstate * (iters + 1) + 1);
            for (int i = 0; i < nstate; i++) {
                wfnat.col(2 * i + 1) = wfna.get().col(i).real(), wfnat.col(2 * i + 2) = wfna.get().col(i).imag();
            } wfnat.col(0) = wfna.getr();
        }

        // calculate the diabatic and adiabatic density matrices and save them
        if (adiabatic) {
            rho = wfna.density(); Pa.row(0)(0) = 0, Pa.row(0).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);
        } rho = wfnd.density(); Pd.row(0)(0) = 0, Pd.row(0).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);

        // print the header
        std::printf("\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIME");

        // perform the propagation
        for (int j = 0; j < iters; j++) {

            // reset the timer, save the previous energy and propagate
            tp = Timer::Now(); double Ep = E; wfnd = wfnd.propagate(R, K);

            // subrtract the lower eigenstates in ITP
            for (int k = 0; k < i && imaginary; k++) wfnd = wfnd - states.at(k) * states.at(k).overlap(wfnd);

            // normalize the wavefunction if ITP and calculate the energy
            if (imaginary) {wfnd = wfnd.normalized();} E = wfnd.energy(Ud);

            // assign the ACF value
            std::complex<double> overlap = states.at(i).overlap(wfnd); acf(j + 1, 0) = (j + 1) * step; acf(j + 1, 1) = overlap.real(); acf(j + 1, 2) = overlap.imag();

            // transform the wfn to the adiabatic basis and save the adiabatic density matrix
            if (adiabatic) {
                wfna = wfnd.adiabatize(UT); rho = wfna.density(); Pa(j + 1, 0) = (j + 1) * step, Pa.row(j + 1).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);
            }

            // save the diabatic density matrix
            rho = wfnd.density(); Pd(j + 1, 0) = (j + 1) * step, Pd.row(j + 1).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);

            // save the diabatic wavefunction to the wavefunction containers
            if (savewfn) {
                for (int k = 0; k < nstate; k++) {
                    wfndt.col(2 * nstate * (j + 1) + 2 * k + 1) = wfnd.get().col(k).real(), wfndt.col(2 * nstate * (j + 1) + 2 * k + 2) = wfnd.get().col(k).imag();
                }
            }

            // save the adiabatic wavefunction to the wavefunction containers
            if (adiabatic && savewfn) {
                for (int k = 0; k < nstate; k++) {
                    wfnat.col(2 * nstate * (j + 1) + 2 * k + 1) = wfna.get().col(k).real(), wfnat.col(2 * nstate * (j + 1) + 2 * k + 2) = wfna.get().col(k).imag();
                }
            }

            // print the iteration info
            std::printf("%6d %20.14f %.2e %s\n", j + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
        }

        // assign results and print a new line
        states.at(i) = wfnd; eps(i) = E; std::cout << std::endl;

        // calculate the spectrum
        Matrix spectrum(iters + 1, 2); spectrum.col(1) = FourierTransform::C2RFFT(acf.col(1) + std::complex<double>(0, 1) * acf.col(2)).array().abs();

        // fill the frequency column of the spectrum
        spectrum.col(0).fill(2 * M_PI / (iters + 1) / step); for (int i = 0; i < iters + 1; i++) spectrum.col(0)(i) *= i - (i < (iters + 1) / 2 ? 0 : (iters + 1));

        // print the populations
        for (int i = 0; i < nstate && adiabatic; i++) {
            std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(nstatessq).reshaped(nstate, nstate).transpose()(i, i));
        } if (adiabatic) std::cout << std::endl;

        // print the populations
        for (int i = 0; i < nstate; i++) {
            std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(nstatessq).reshaped(nstate, nstate).transpose()(i, i));
        } std::cout << std::endl;

        // save the resulting data data
        MEASURE("WAVEFUNCTIONS, DENSITY MATRICES, ACF, SPECTRUM AND POTENTIAL WRITING: ",
            if (             savewfn) Eigen::Write("PSI_DIA_"      + std::to_string(i) + ".mat", wfndt   );
                                      Eigen::Write("P_DIA_"        + std::to_string(i) + ".mat", Pd      );
            if (adiabatic           ) Eigen::Write("U_ADIA_"       + std::to_string(i) + ".mat", Ua      );
            if (adiabatic && savewfn) Eigen::Write("PSI_ADIA_"     + std::to_string(i) + ".mat", wfnat   );
            if (adiabatic           ) Eigen::Write("P_ADIA_"       + std::to_string(i) + ".mat", Pa      );
                                      Eigen::Write("ACF_DIA_"      + std::to_string(i) + ".mat", acf     );
                                      Eigen::Write("SPECTRUM_DIA_" + std::to_string(i) + ".mat", spectrum);
        )
    }

    // print the state energies
    std::cout << std::endl; for (int i = 0; i < program.get<int>("-o"); i++) std::printf("EIGENSTATE %d ENERGY: %.14f\n", i, eps(i));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

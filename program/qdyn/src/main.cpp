#include "fourier.h"
#include "timer.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Quantum Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-d", "--dimension").help("-- Dimension of the provided potential and wavefunction.").default_value(1).scan<'i', int>();
    program.add_argument("-f", "--factor").help("-- Factor for scaling the saved wavefunction.").default_value(1.0).scan<'g', double>();
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-o", "--optimize").help("-- Perform the dynamics in imaginary time for several orthogonal states.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("--adiabatic").help("-- Enable transform to adiabatic basis.").default_value(false).implicit_value(true);
    program.add_argument("--align").help("-- Align the wavefunction values to the potential.").default_value(false).implicit_value(true);
    program.add_argument("--savewfn").help("-- Save the time evolution of the wavefunction.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    int dim = program.get<int>("-d"), iters = program.get<int>("-i"), optstates = program.get<int>("-o"); double factor = program.get<double>("-f"), mass = program.get<double>("-m"), step = program.get<double>("-s"), momentum = program.get<double>("-p");

    // extract boolean flags
    bool adiabatic = program.get<bool>("--adiabatic"), align = program.get<bool>("--align"), imaginary = program.is_used("--optimize"), savewfn = program.get<bool>("--savewfn");

    // load the potential and initial wavefunction
    MEASURE("READING POTENTIAL AND INITIAL WAVEFUNCTION: ", Matrix Ud = Eigen::LoadMatrix("U_DIA.mat"); 
        std::vector<Wavefunction> states(optstates, Wavefunction(Eigen::LoadMatrix("PSI_DIA_GUESS.mat").rightCols(2 * (int)std::sqrt(Ud.cols() - dim)), Ud.leftCols(dim), mass, momentum).normalized());
    )

    // extract the number of states
    int nstate = states.at(0).get().cols(), nstatessq = nstate * nstate;

    // remove the first column from matrix (the independent variable column)
    Ud.block(0, 0, Ud.rows(), Ud.cols() - dim) = Ud.rightCols(Ud.cols() - dim); Ud.conservativeResize(Ud.rows(), Ud.cols() - dim);

    // define and initialize the adiabatic variables
    Matrix Ua; std::vector<Matrix> UT; if (adiabatic) {
        Ua = Matrix(Ud.rows(), nstate + dim); UT = std::vector<Matrix>(Ud.rows(), Matrix::Identity(nstate, nstate));
    }

    // diagonalize the potential at each point
    for (int i = 0; i < Ud.rows() && adiabatic; i++) {

        // initialize the diabatic potential at the current point
        Matrix UD = Ud.row(i).reshaped(nstate, nstate).transpose();

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Matrix> solver(UD); Matrix C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(nstate, 1);

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
    if (adiabatic) {MEASURE("\nSAVING THE ADIABATIC POTENTIAL: ", Eigen::Write("U_ADIA.mat", Ua))}

    // define the real time kinetic and potential operators in diabatic basis
    auto[R, K] = states.at(0).propagator(Ud, std::complex<double>(imaginary, !imaginary), step);

    // define the energy, wfn and acf vector
    Vector eps(optstates); Matrix acf(imaginary ? 0 : iters + 1, imaginary ? 0 : 3);

    // for every orthogonal state
    for (int i = 0; i < optstates; i++) {

        // set the initial wfn, energy and acf value
        Wavefunction wfnd = states.at(i); double E = wfnd.energy(Ud); if (!imaginary) acf(0, 0) = 0, acf(0, 1) = 1, acf(0, 2) = 0;

        // define the diabatic and adiabatic wavefunction and density matrix containers
        Wavefunction wfna; Matrix wfndt, wfnat, Pd, Pa, rho; if (nstate > 1) Pd = Matrix(iters + 1, nstatessq + 1), Pa = Matrix(iters + 1, nstatessq + 1), rho = Matrix(nstate, nstate);

        // transform the wfn to the adiabatic basis
        if (adiabatic) wfna = wfnd.adiabatize(UT);

        // save the diabatic wavefunction and fill the independent variable columns
        if (savewfn) {
            wfndt = Matrix(Ud.rows(), 2 * nstate * (iters + 1) + dim);
            for (int i = 0; i < nstate; i++) {
                wfndt.col(2 * i + dim) = factor * wfnd.get().col(i).real(), wfndt.col(2 * i + dim + 1) = factor * wfnd.get().col(i).imag();
                if (align) {
                    wfndt.col(2 * i + dim + 0) = wfndt.col(2 * i + dim + 0).array() + Ud.col(i * (nstate + 1)).array();
                    wfndt.col(2 * i + dim + 1) = wfndt.col(2 * i + dim + 1).array() + Ud.col(i * (nstate + 1)).array();
                }
            } wfndt.leftCols(dim) = wfnd.getr();
        }

        // save the adiabatic wavefunction and fill the independent variable columns
        if (adiabatic && savewfn) {
            wfnat = Matrix(Ua.rows(), 2 * nstate * (iters + 1) + dim);
            for (int i = 0; i < nstate; i++) {
                wfnat.col(2 * i + dim) = factor * wfna.get().col(i).real(), wfnat.col(2 * i + dim + 1) = factor * wfna.get().col(i).imag();
                if (align) {
                    wfnat.col(2 * i + dim + 0) = wfnat.col(2 * i + dim + 0).array() + Ua.col(1 + i).array();
                    wfnat.col(2 * i + dim + 1) = wfnat.col(2 * i + dim + 1).array() + Ua.col(1 + i).array();
                }
            } wfnat.leftCols(dim) = wfna.getr();
        }

        // calculate the diabatic and adiabatic density matrices and save them
        if (nstate > 1) {
            if (adiabatic) {
                rho = wfna.density(); Pa.row(0)(0) = 0, Pa.row(0).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);
            } rho = wfnd.density(); Pd.row(0)(0) = 0, Pd.row(0).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);
        }

        // print the header
        std::printf("\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIME");

        // perform the propagation
        for (int j = 0; j < iters; j++) {

            // reset the timer, save the previous energy and propagate
            Timer::Timepoint tsit = Timer::Now(); double Ep = E; wfnd = wfnd.propagate(R, K);

            // subrtract the lower eigenstates in ITP
            for (int k = 0; k < i && imaginary; k++) wfnd = wfnd - states.at(k) * states.at(k).overlap(wfnd);

            // normalize the wavefunction if ITP and calculate the energy
            if (imaginary) {wfnd = wfnd.normalized();} E = wfnd.energy(Ud);

            // assign the ACF value
            std::complex<double> overlap; if (!imaginary) {overlap = states.at(i).overlap(wfnd); acf(j + 1, 0) = (j + 1) * step; acf(j + 1, 1) = overlap.real(); acf(j + 1, 2) = overlap.imag();}

            // transform the wfn to the adiabatic basis and save the adiabatic density matrix
            if (adiabatic) {
                wfna = wfnd.adiabatize(UT); if (nstate > 1) {rho = wfna.density(); Pa(j + 1, 0) = (j + 1) * step, Pa.row(j + 1).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);}
            }

            // save the diabatic density matrix
            if (nstate > 1) {rho = wfnd.density(); Pd(j + 1, 0) = (j + 1) * step, Pd.row(j + 1).rightCols(nstatessq) = rho.transpose().reshaped(1, nstatessq);}

            // save the diabatic wavefunction to the wavefunction containers
            if (savewfn) {
                for (int k = 0; k < nstate; k++) {
                    wfndt.col(2 * nstate * (j + 1) + 2 * k + dim) = factor * wfnd.get().col(k).real(), wfndt.col(2 * nstate * (j + 1) + 2 * k + dim + 1) = factor * wfnd.get().col(k).imag();
                    if (align) {
                        wfndt.col(2 * nstate * (j + 1) + 2 * k + dim + 0) = wfndt.col(2 * nstate * (j + 1) + 2 * k + dim + 0).array() + Ud.col(k * (nstate + 1)).array();
                        wfndt.col(2 * nstate * (j + 1) + 2 * k + dim + 1) = wfndt.col(2 * nstate * (j + 1) + 2 * k + dim + 1).array() + Ud.col(k * (nstate + 1)).array();
                    }
                }
            }

            // save the adiabatic wavefunction to the wavefunction containers
            if (adiabatic && savewfn) {
                for (int k = 0; k < nstate; k++) {
                    wfnat.col(2 * nstate * (j + 1) + 2 * k + dim) = factor * wfna.get().col(k).real(), wfnat.col(2 * nstate * (j + 1) + 2 * k + dim + 1) = factor * wfna.get().col(k).imag();
                    if (align) {
                        wfnat.col(2 * nstate * (j + 1) + 2 * k + dim + 0) = wfnat.col(2 * nstate * (j + 1) + 2 * k + dim + 0).array() + Ua.col(1 + k).array();
                        wfnat.col(2 * nstate * (j + 1) + 2 * k + dim + 1) = wfnat.col(2 * nstate * (j + 1) + 2 * k + dim + 1).array() + Ua.col(1 + k).array();
                    }
                }
            }

            // print the iteration info
            std::printf("%6d %20.14f %.2e %s\n", j + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tsit)).c_str());
        }

        // assign results and print a new line
        states.at(i) = wfnd; eps(i) = E; std::cout << std::endl;

        // define the spectrum matrix
        Matrix spectrum(imaginary ? 0 : iters + 1, imaginary ? 0 : 2);

        // calculate the spectrum
        if (!imaginary) {

            // perform the fourier transform of the ACF
            spectrum.col(1) = FourierTransform::C2RFFT(acf.col(1) + std::complex<double>(0, 1) * acf.col(2)).array().abs();

            // fill the frequency column of the spectrum
            spectrum.col(0).fill(2 * M_PI / (iters + 1) / step); for (int i = 0; i < iters + 1; i++) spectrum.col(0)(i) *= i - (i < (iters + 1) / 2 ? 0 : (iters + 1));
        }

        // print the populations
        for (int i = 0; i < nstate && nstate > 1; i++) {
            std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(nstatessq).reshaped(nstate, nstate).transpose()(i, i));
        } if (nstate > 1) std::cout << std::endl;

        // print the populations
        for (int i = 0; i < nstate && nstate > 1 && adiabatic; i++) {
            std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(nstatessq).reshaped(nstate, nstate).transpose()(i, i));
        } if (nstate > 1 && adiabatic) std::cout << std::endl;

        // save the resulting data data
        MEASURE("WAVEFUNCTIONS, DENSITY MATRICES, ACF AND SPECTRUM WRITING: ",
            if (             savewfn) Eigen::Write("PSI_DIA_"      + std::to_string(i) + ".mat", wfndt   );
            if (nstate > 1          ) Eigen::Write("P_DIA_"        + std::to_string(i) + ".mat", Pd      );
            if (adiabatic && savewfn) Eigen::Write("PSI_ADIA_"     + std::to_string(i) + ".mat", wfnat   );
            if (adiabatic           ) Eigen::Write("P_ADIA_"       + std::to_string(i) + ".mat", Pa      );
            if (!imaginary          ) Eigen::Write("ACF_DIA_"      + std::to_string(i) + ".mat", acf     );
            if (!imaginary          ) Eigen::Write("SPECTRUM_DIA_" + std::to_string(i) + ".mat", spectrum);
        )
    }

    // print the state energies
    std::cout << std::endl; for (int i = 0; i < optstates; i++) std::printf("FINAL WFN %02d ENERGY: %.14f\n", i, eps(i));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

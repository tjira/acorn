#include "timer.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Nonadiabatic Quantum Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
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
    bool adiabatic = program.get<bool>("--adiabatic"), savewfn = program.get<bool>("--savewfn");

    // load the potential in the diabatic basis
    MEASURE("INITIAL WAVEFUNCTION AND POTENTIAL IN DIABATIC BASIS READING: ",
        Matrix Ud = Eigen::LoadMatrix("U_DIA.mat"); Wavefunction wfnd(Eigen::LoadMatrix("PSI_DIA_GUESS.mat").rightCols(2 * (int)std::sqrt(Ud.cols() - 1)), Ud.leftCols(1), mass, momentum);
    )

    // extract the number of states
    int states = wfnd.states(), statessq = states * states;

    // normalize the wfn, remove first column from matrix (the independent variable column)
    wfnd = wfnd.normalized(); Ud.block(0, 0, Ud.rows(), Ud.cols() - 1) = Ud.rightCols(Ud.cols() - 1); Ud.conservativeResize(Ud.rows(), Ud.cols() - 1);

    // define wavefunction value containers, density matrices and adiabatic transformation matrices
    Matrix wfndt, wfnat, Pd(iters + 1, statessq + 1), Pa, Ua(Ud.rows(), states + 1), rho(states, states); std::vector<Matrix> UT; if (savewfn) wfndt = Matrix(Ud.rows(), 2 * states * (iters + 1) + 1);

    // initialize the adiabatic variables
    if (adiabatic) {
        Ua = Matrix(Ud.rows(), states + 1), Pa = Matrix(iters + 1, statessq + 1);
        UT = std::vector<Matrix>(Ud.rows(), Matrix::Identity(states, states));
        if (savewfn) {wfnat = Matrix(wfndt.rows(), wfndt.cols());}
    }

    // diagonalize the potential at each point
    for (int i = 0; i < Ud.rows() && adiabatic; i++) {

        // initialize the diabatic potential at the current point
        Matrix UD = Ud.row(i).reshaped(states, states).transpose();

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Matrix> solver(UD); Matrix C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(states, 1);

        // fill the diabatic matrix
        Ua(i, 0) = wfnd.getr()(i); Ua.row(i).rightCols(states) = eps.transpose();

        // calculate the overlap of eigenvectors
        for (int j = 0; j < states; j++) {
            overlap(j) = UT.at(std::max(i - 1, 0)).col(j).transpose() * C.col(j); overlap(j) = overlap(j) < 0 ? -1 : 1;
        }

        // maximize the overlap with the previuos transformation
        UT.at(i) = C * overlap.asDiagonal();
    }

    // adiabatize the guess wavefunction and fill the independent variable columns
    Wavefunction wfna; if (adiabatic) {wfna = wfnd.adiabatize(UT); if (savewfn) wfnat.col(0) = wfna.getr();} if (savewfn) wfndt.col(0) = wfnd.getr();

    // calculate the diabatic and adiabatic density matrices and save them
    if (adiabatic) {
        rho = wfna.density(); Pa.row(0)(0) = 0, Pa.row(0).rightCols(statessq) = rho.transpose().reshaped(1, statessq);
    } rho = wfnd.density(); Pd.row(0)(0) = 0, Pd.row(0).rightCols(statessq) = rho.transpose().reshaped(1, statessq);

    // fill the initial wavefunction values of the adiabatic evolution matrices
    if (adiabatic && savewfn) {
        for (int i = 0; i < states; i++) {
            wfnat.col(2 * i + 1) = wfna.get().col(i).real(), wfnat.col(2 * i + 2) = wfna.get().col(i).imag();
        }
    }

    std::cout << states << std::endl;

    // fill the initial wavefunction values of the diabatic evolution matrices
    if (savewfn) {
        for (int i = 0; i < states; i++) {
            wfndt.col(2 * i + 1) = wfnd.get().col(i).real(), wfndt.col(2 * i + 2) = wfnd.get().col(i).imag();
        }
    }

    // define the real time kinetic and potential operators in diabatic basis and the energy container
    auto[R, K] = wfnd.propagator(Ud, std::complex<double>(0, 1), step); double E = wfnd.energy(Ud);

    // print the header
    std::printf("\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIME");

    // perform the propagation
    for (int i = 0; i < iters; i++) {

        // reset the timer and save the previous energy
        tp = Timer::Now(); double Ep = E;

        // apply the propagator, calculate energy and transform to the adiabatic basis
        wfnd = wfnd.propagate(R, K), E = wfnd.energy(Ud); if (adiabatic) wfna = wfnd.adiabatize(UT);

        // save the adiabatic density matrix
        if (adiabatic) {
            rho = wfna.density(); Pa(i + 1, 0) = (i + 1) * step, Pa.row(i + 1).rightCols(statessq) = rho.transpose().reshaped(1, statessq);
        }

        // save the diabatic density matrix
        rho = wfnd.density(); Pd(i + 1, 0) = (i + 1) * step, Pd.row(i + 1).rightCols(statessq) = rho.transpose().reshaped(1, statessq);

        // save the adiabatic wavefunction to the wavefunction containers
        if (adiabatic && savewfn) {
            for (int j = 0; j < states; j++) {
                wfnat.col(2 * states * (i + 1) + 2 * j + 1) = wfna.get().col(j).real(), wfnat.col(2 * states * (i + 1) + 2 * j + 2) = wfna.get().col(j).imag();
            }
        }

        // save the diabatic wavefunction to the wavefunction containers
        if (savewfn) {
            for (int j = 0; j < states; j++) {
                wfndt.col(2 * states * (i + 1) + 2 * j + 1) = wfnd.get().col(j).real(), wfndt.col(2 * states * (i + 1) + 2 * j + 2) = wfnd.get().col(j).imag();
            }
        }

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
    }

    // print a new line
    std::cout << std::endl;

    // print the populations
    for (int i = 0; i < states && adiabatic; i++) {
        std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(statessq).reshaped(states, states).transpose()(i, i));
    } std::cout << std::endl;

    // print the populations
    for (int i = 0; i < states; i++) {
        std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(statessq).reshaped(states, states).transpose()(i, i));
    } std::cout << std::endl;

    // save the resulting data data
    MEASURE("WAVEFUNCTIONS, DENSITY MATRICES AND POTENTIAL WRITING: ",
        if (             savewfn) Eigen::Write("PSI_DIA.mat" , wfndt);
                                  Eigen::Write("P_DIA.mat"   , Pd   );
        if (adiabatic           ) Eigen::Write("U_ADIA.mat"  , Ua   );
        if (adiabatic && savewfn) Eigen::Write("PSI_ADIA.mat", wfnat);
        if (adiabatic           ) Eigen::Write("P_ADIA.mat"  , Pa   );
    )

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

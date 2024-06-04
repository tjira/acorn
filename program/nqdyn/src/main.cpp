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
        Matrix Ud = Eigen::LoadMatrix("U_DIA.mat"); Wavefunction<2> wfnd(Eigen::LoadMatrix("PSI_DIA_GUESS.mat").rightCols(4), Ud.leftCols(1), mass, momentum);
    )

    // normalize the wfn, remove first column from matrix (the independent variable column)
    wfnd = wfnd.normalized(); Ud.block(0, 0, Ud.rows(), Ud.cols() - 1) = Ud.rightCols(Ud.cols() - 1); Ud.conservativeResize(Ud.rows(), Ud.cols() - 1);

    // define wavefunction value containers, density matrices and adiabatic transformation matrices
    Matrix wfndt, wfnat, Pd(iters + 1, 5), Pa, Ua(Ud.rows(), 3), rho(2, 2); std::vector<Matrix> UT; if (savewfn) wfndt = Matrix(Ud.rows(), 5 + 4 * iters);

    // initialize the adiabatic variables
    if (adiabatic) {
        UT = std::vector<Matrix>(Ud.rows(), Matrix::Identity(2, 2));
        Ua = Matrix(Ud.rows(), 3), Pa = Matrix(iters + 1, 5);
        if (savewfn) {wfnat = Matrix(Ud.rows(), 5 + 4 * iters);}
    }

    // diagonalize the potential at each point
    for (int i = 0; i < Ud.rows() && adiabatic; i++) {

        // define the diabatic potential at the current point and fill it
        Matrix UD(2, 2); UD << Ud(i, 0), Ud(i, 1), Ud(i, 2), Ud(i, 3);

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Matrix> solver(UD); Matrix C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(2, 1);

        // fill the diabatic matrix
        Ua.row(i) << wfnd.getr()(i), eps(0), eps(1);

        // calculate the overlap of eigenvectors
        for (size_t j = 0; j < 2; j++) {
            overlap(j) = UT.at(std::max(i - 1, 0)).col(j).transpose() * C.col(j); overlap(j) /= std::abs(overlap(j));
        }

        // maximize the overlap with the previuos transformation
        UT.at(i) = C * overlap.asDiagonal();
    }

    // adiabatize the guess wavefunction and fill the independent variable columns
    Wavefunction<2> wfna; if (adiabatic) {wfna = wfnd.adiabatize(UT); if (savewfn) wfnat.col(0) = wfna.getr();} if (savewfn) wfndt.col(0) = wfnd.getr();

    // calculate the diabatic and adiabatic density matrices and save them
    if (adiabatic) {
        rho = wfna.density(); Pa.row(0) << 0, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);
    } rho = wfnd.density(); Pd.row(0) << 0, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);

    // fill the initial wavefunction values of the adiabatic evolution matrices
    if (adiabatic && savewfn) {
        wfnat.col(1) = wfna.get().col(0).real(), wfnat.col(2) = wfna.get().col(0).imag();
        wfnat.col(3) = wfna.get().col(1).real(), wfnat.col(4) = wfna.get().col(1).imag();
    }

    // fill the initial wavefunction values of the diabatic evolution matrices
    if (savewfn) wfndt.col(1) = wfnd.get().col(0).real(), wfndt.col(2) = wfnd.get().col(0).imag();
    if (savewfn) wfndt.col(3) = wfnd.get().col(1).real(), wfndt.col(4) = wfnd.get().col(1).imag();

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
            rho = wfna.density(); Pa.row(i + 1) << (i + 1) * step, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);
        }

        // save the diabatic density matrix
        rho = wfnd.density(); Pd.row(i + 1) << (i + 1) * step, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);

        // save the adiabatic wavefunction to the wavefunction containers
        if (adiabatic && savewfn) {
            wfnat.col(5 + 4 * i) = wfna.get().col(0).real(), wfnat.col(6 + 4 * i) = wfna.get().col(0).imag();
            wfnat.col(7 + 4 * i) = wfna.get().col(1).real(), wfnat.col(8 + 4 * i) = wfna.get().col(1).imag();
        }

        // save the diabatic wavefunction to the wavefunction containers
        if (savewfn) wfndt.col(5 + 4 * i) = wfnd.get().col(0).real(), wfndt.col(6 + 4 * i) = wfnd.get().col(0).imag();
        if (savewfn) wfndt.col(7 + 4 * i) = wfnd.get().col(1).real(), wfndt.col(8 + 4 * i) = wfnd.get().col(1).imag();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
    }

    // print a new line
    std::cout << std::endl;

    // print the populations
    for (int i = 0; i < (int)std::sqrt(Pd.cols() - 1); i++) {
        std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(4).reshaped(2, 2)(i, i));
    } std::cout << std::endl;

    // print the populations
    for (int i = 0; i < (int)std::sqrt(Pd.cols() - 1); i++) {
        std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(4).reshaped(2, 2)(i, i));
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

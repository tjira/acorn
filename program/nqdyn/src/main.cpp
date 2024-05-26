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

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters
    int iters = program.get<int>("-i"); double mass = program.get<double>("-m"), step = program.get<double>("-s"), momentum = program.get<double>("-p");

    // load the potential in the diabatic basis
    MEASURE("DIABATIC POTENTIAL MATRIX READING: ", EigenMatrix<> Ud = Eigen::LoadMatrix("U_DIA.mat"))

    // load the wavefunction in the diabatic basis
    MEASURE("DIABATIC WAVEFUNCTION READING:     ", Wavefunction<2> wfnd(Eigen::LoadMatrix("PSI_DIA_GUESS.mat").rightCols(4), Ud.leftCols(1), mass, momentum);)

    // normalize the wfn, remove first column from matrix (the independent variable column)
    wfnd = wfnd.normalized(); Ud.block(0, 0, Ud.rows(), Ud.cols() - 1) = Ud.rightCols(Ud.cols() - 1); Ud.conservativeResize(Ud.rows(), Ud.cols() - 1);

    // define wavefunction value containers, density matrices and adiabatic transformation matrices
    EigenMatrix<> wfndt(Ud.rows(), 5 + 4 * iters), wfnat(Ud.rows(), 5 + 4 * iters), Pd(iters + 1, 5), Pa(iters + 1, 5);
    std::vector<EigenMatrix<>> UT(Ud.rows(), EigenMatrix<>::Identity(2, 2)); EigenMatrix<> Ua(Ud.rows(), 3), rho(2, 2);

    // diagonalize the potential at each point
    for (int i = 0; i < Ud.rows(); i++) {

        // define the diabatic potential at the current point and fill it
        EigenMatrix<> UD(2, 2); UD << Ud(i, 0), Ud(i, 1), Ud(i, 2), Ud(i, 3);

        // solve the eigenvalue problem and define overlap
        auto[eps, C] = Eigen::Saes(UD); EigenVector<> overlap(2);

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
    Wavefunction<2> wfna = wfnd.adiabatize(UT); wfndt.col(0) = wfnd.getr(), wfnat.col(0) = wfna.getr();

    // calculate the diabatic and adiabatic density matrices and save them
    rho = wfnd.density(); Pd.row(0) << 0, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);
    rho = wfna.density(); Pa.row(0) << 0, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);

    // fill the initial wavefunction values of the evolution matrices
    wfnat.col(1) = wfna.get().col(0).real(), wfnat.col(2) = wfna.get().col(0).imag();
    wfnat.col(3) = wfna.get().col(1).real(), wfnat.col(4) = wfna.get().col(1).imag();
    wfndt.col(1) = wfnd.get().col(0).real(), wfndt.col(2) = wfnd.get().col(0).imag();
    wfndt.col(3) = wfnd.get().col(1).real(), wfndt.col(4) = wfnd.get().col(1).imag();

    // define the real time kinetic and potential operators in diabatic basis and the energy container
    auto[R, K] = wfnd.propagator(Ud, std::complex<double>(0, 1), step); double E = wfnd.energy(Ud);

    // print the header
    std::printf("\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIME");

    // perform the propagation
    for (int i = 0; i < iters; i++) {

        // reset the timer and save the previous energy
        tp = Timer::Now(); double Ep = E;

        // apply the propagator, calculate energy and transform to the adiabatic basis
        wfnd = wfnd.propagate(R, K), E = wfnd.energy(Ud); wfna = wfnd.adiabatize(UT);

        rho = wfnd.density(); Pd.row(i + 1) << (i + 1) * step, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);
        rho = wfna.density(); Pa.row(i + 1) << (i + 1) * step, rho(0, 0), rho(0, 1), rho(1, 0), rho(1, 1);

        // save the diabatic and adiabatic wavefunction to the wavefunction containers
        wfnat.col(5 + 4 * i) = wfna.get().col(0).real(), wfnat.col(6 + 4 * i) = wfna.get().col(0).imag();
        wfnat.col(7 + 4 * i) = wfna.get().col(1).real(), wfnat.col(8 + 4 * i) = wfna.get().col(1).imag();
        wfndt.col(5 + 4 * i) = wfnd.get().col(0).real(), wfndt.col(6 + 4 * i) = wfnd.get().col(0).imag();
        wfndt.col(7 + 4 * i) = wfnd.get().col(1).real(), wfndt.col(8 + 4 * i) = wfnd.get().col(1).imag();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
    }

    // print a new line
    std::cout << std::endl;

    // save the diabatic data
    MEASURE("DIABATIC DENSITY MATRIX WRITING: ", Eigen::Write("P_DIA.mat", Pd))
    MEASURE("DIABATIC WAVEFUNCTION WRITING:   ", Eigen::Write("PSI_DIA.mat", wfndt))

    // print a new line
    std::cout << std::endl;

    // save the adiabatic data
    MEASURE("ADIABATIC DENSITY MATRIX WRITING: ", Eigen::Write("P_ADIA.mat", Pa))
    MEASURE("ADIABATIC WAVEFUNCTION WRITING:   ", Eigen::Write("PSI_ADIA.mat", wfnat))
    MEASURE("ADIABATIC POTENTIAL WRITING:      ", Eigen::Write("U_ADIA.mat", Ua))

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

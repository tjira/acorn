#include "timer.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Nonadiabatic Quantum Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(350).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(2000.0).scan<'g', double>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(10.95).scan<'g', double>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(10.0).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters
    int iters = program.get<int>("-i"); double mass = program.get<double>("-m"), step = program.get<double>("-s");

    // extract the guess expressions to an array
    std::array<std::string, 2> guess = {program.get<std::vector<std::string>>("-g").at(0), program.get<std::vector<std::string>>("-g").at(1)};

    // load the potential
    tp = Timer::Now(); std::cout << "POTENTIAL READING: " << std::flush; EigenMatrix<> U = Eigen::LoadMatrix("U.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // load the wavefunction
    tp = Timer::Now(); std::cout << "WFUNCTION READING: " << std::flush; Wavefunction<2> wfn(Eigen::LoadMatrix("PSI_0.mat").rightCols(4), U.leftCols(1), mass, program.get<double>("-p")); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // normalize the wfn, remove first column from matrix (the independent variable column)
    wfn = wfn.normalized(); U.block(0, 0, U.rows(), U.cols() - 1) = U.rightCols(U.cols() - 1); U.conservativeResize(U.rows(), U.cols() - 1);

    // define the density matrix and calculate its initial values
    EigenMatrix<> density = wfn.density(); EigenMatrix<> P(iters + 1, 5); P.row(0) << 0.0, density(0, 0), density(0, 1), density(1, 0), density(1, 1);

    // define the real time kinetic and potential operators, energy and density matrix container
    auto[R, K] = wfn.propagator(U, std::complex<double>(0, 1), step); double E = wfn.energy(U);

    // print the header
    std::printf("\n%6s %20s %20s %20s %8s %12s\n", "ITER", "ENERGY", "D00", "D11", "|dE|", "TIME");

    // perform the propagation
    for (int i = 0; i < iters; i++) {

        // reset the timer, save the previous values of wfn and energy
        tp = Timer::Now(); Wavefunction<2> wfnp = wfn; double Ep = E;

        // apply the propagator
        wfn = wfn.propagate(R, K);

        // calculate the total energy and density matrix
        E = wfn.energy(U), density = wfn.density(); P.row(i + 1) << ((i + i) * step), density(0, 0), density(0, 1), density(1, 0), density(1, 1);

        // print the iteration info
        std::printf("%6d %20.14f %20.14f %20.14f %.2e %s\n", i + 1, E, density(0, 0), density(1, 1), std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
    }

    // print new line
    std::cout << std::endl;

    // save the density matrix to a file
    tp = Timer::Now(); std::cout << "P MATRIX WRITING: " << std::flush; Eigen::Write("P.mat", P); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

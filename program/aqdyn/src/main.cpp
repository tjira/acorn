#include "timer.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Adiabatic Quantum Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(1000).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-n", "--nstate").help("-- Number of states to calculate.").default_value(1).scan<'i', int>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(0.1).scan<'g', double>();
    program.add_argument("--imaginary").help("-- Enable the imaginary time propagation.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters
    int iters = program.get<int>("-i"), nstate = program.get<int>("-n"); double mass = program.get<double>("-m"), step = program.get<double>("-s");

    // load the potential
    tp = Timer::Now(); std::cout << "POTENTIAL READING: " << std::flush; EigenMatrix<> U = Eigen::LoadMatrix("U.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // load the wavefunction
    tp = Timer::Now(); std::cout << "WFUNCTION READING: " << std::flush; Wavefunction<1> wfn(Eigen::LoadMatrix("PSI0_0.mat").rightCols(2), U.leftCols(U.cols() - 1), mass); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // normalize the wfn, extract the potential values from the last column
    wfn = wfn.normalized(); U.col(0) = U.rightCols(1); U.conservativeResize(U.rows(), 1);

    // calculate the propagator
    auto[R, K] = wfn.propagator(U, std::complex<double>(program.get<bool>("--imaginary"), !program.get<bool>("--imaginary")), step);

    // define the energy and wfn vector
    EigenVector<> eps(nstate); std::vector<Wavefunction<1>> states(nstate, wfn);

    // perform the dynamics for every state
    for (int i = 0; i < nstate; i++) {

        // assign the initial wfn with its energy and create the evolution matrix
        wfn = states.at(i); double E = wfn.energy(U); EigenMatrix<> evolution(U.rows(), 3 + 2 * iters);

        // save the initial state
        evolution.col(0) = wfn.getr(); evolution.col(1) = wfn.get().at(0).real(), evolution.col(2) = wfn.get().at(0).imag();

        // print the header
        std::printf("\nSTATE %d\n%6s %20s %8s %12s\n", i, "ITER", "ENERGY", "|dE|", "TIME");

        // propagate the current state
        for (int j = 0; j < iters; j++) {

            // reset the timer, save the energy and propagate
            tp = Timer::Now(); double Ep = E; wfn = wfn.propagate(R, K);

            // subrtract the lower eigenstates
            for (int k = 0; k < i && program.get<bool>("--imaginary"); k++) wfn = wfn - states.at(k) * wfn.overlap(states.at(k));

            // normalize the wavefunction and calculate the energy
            wfn = wfn.normalized(); E = wfn.energy(U);

            // save the wavefunction to the evolution matrix
            evolution.col(3 + 2 * j) = wfn.get().at(0).real(); evolution.col(4 + 2 * j) = wfn.get().at(0).imag();

            // print the iteration info
            std::printf("%6d %20.14f %.2e %s\n", j + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
        }

        // append to states and energy, print the new line and save the evolution
        states.at(i) = wfn; eps(i) = E; Eigen::Write("PSI" + std::to_string(i) + "_T.mat", evolution);
    }

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

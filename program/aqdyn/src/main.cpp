#include "fourier.h"
#include "timer.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Adiabatic Quantum Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-n", "--nstate").help("-- Number of states to calculate.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("--imaginary").help("-- Enable the imaginary time propagation.").default_value(false).implicit_value(true);
    program.add_argument("--savewfn").help("-- Save the time evolution of the wavefunction.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters
    int iters = program.get<int>("-i"), nstate = program.get<int>("-n"); double mass = program.get<double>("-m"), step = program.get<double>("-s");

    // extract boolean flags
    bool imaginary = program.get<bool>("--imaginary"), savewfn = program.get<bool>("--savewfn");

    // load the initial wavefunction and potential
    MEASURE("INITIAL WAVEFUNCTION AND POTENTIAL READING: ",
        Matrix U = Eigen::LoadMatrix("U_ADIA.mat"); Wavefunction<1> wfn(Eigen::LoadMatrix("PSI_ADIA_GUESS.mat").rightCols(2), U.leftCols(U.cols() - 1), mass, program.get<double>("-p"));
    )

    // normalize the wfn, extract the potential values from the last column
    wfn = wfn.normalized(); U.col(0) = U.rightCols(1); U.conservativeResize(U.rows(), 1);

    // calculate the propagator
    auto[R, K] = wfn.propagator(U, std::complex<double>(imaginary, !imaginary), step);

    // define the energy, wfn and acf vector
    Vector eps(nstate); Matrix acf(iters + 1, 3); std::vector<Wavefunction<1>> states(nstate, wfn);

    // perform the dynamics for every state
    for (int i = 0; i < nstate; i++) {

        // assign the initial wfn with its energy and acf value
        wfn = states.at(i); double E = wfn.energy(U); acf(0, 0) = 0, acf(0, 1) = 1, acf(0, 2) = 0;

        // save the initial state
        Matrix wfnt; if (savewfn) {
            wfnt = Matrix(U.rows(), 3 + 2 * iters); wfnt.col(0) = wfn.getr(); wfnt.col(1) = wfn.get().col(0).real(), wfnt.col(2) = wfn.get().col(0).imag();
        }

        // print the header
        std::printf("\nSTATE %d\n%6s %20s %8s %12s\n", i, "ITER", "ENERGY", "|dE|", "TIME");

        // propagate the current state
        for (int j = 0; j < iters; j++) {

            // reset the timer, save the energy and propagate
            tp = Timer::Now(); double Ep = E; wfn = wfn.propagate(R, K);

            // subrtract the lower eigenstates
            for (int k = 0; k < i && imaginary; k++) wfn = wfn - states.at(k) * wfn.overlap(states.at(k));

            // normalize the wavefunction and calculate the energy
            if (imaginary) {wfn = wfn.normalized();} E = wfn.energy(U);

            // assign the acf value
            std::complex<double> overlap = wfn.overlap(states.at(i)); acf(j + 1, 0) = (j + 1) * step; acf(j + 1, 1) = overlap.real(); acf(j + 1, 2) = overlap.imag();

            // save the wavefunction to the evolution matrix
            if (savewfn) {wfnt.col(3 + 2 * j) = wfn.get().col(0).real(); wfnt.col(4 + 2 * j) = wfn.get().col(0).imag();}

            // print the iteration info
            std::printf("%6d %20.14f %.2e %s\n", j + 1, E, std::abs(E - Ep), Timer::Format(Timer::Elapsed(tp)).c_str());
        }

        // print the new line
        std::cout << std::endl;

        // append to states and energy, print the new line and save the evolution
        states.at(i) = wfn; eps(i) = E; if (savewfn) {MEASURE("SAVING THE WFN OF STATE " + std::to_string(i) + ":      ", Eigen::Write("PSI_ADIA_" + std::to_string(i) + ".mat", wfnt))}

        // save the acf
        MEASURE("SAVING THE ACF OF STATE " + std::to_string(i) + ":      ", Eigen::Write("ACF_ADIA_" + std::to_string(i) + ".mat", acf))

        // calculate the spectrum
        Matrix spectrum(iters + 1, 2); spectrum.col(1) = FourierTransform::C2RFFT(acf.col(1) + std::complex<double>(0, 1) * acf.col(2)).array().abs();

        // fill the frequency column
        spectrum.col(0).fill(2 * M_PI / (iters + 1) / step); for (int i = 0; i < iters + 1; i++) spectrum.col(0)(i) *= i - (i < (iters + 1) / 2 ? 0 : (iters + 1));

        // save the spectrum
        MEASURE("SAVING THE SPECTRUM OF STATE " + std::to_string(i) + ": ", Eigen::Write("SPECTRUM_ADIA_" + std::to_string(i) + ".mat", spectrum))
    }

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

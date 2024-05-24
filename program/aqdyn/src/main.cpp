#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Program", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program.add_argument("-g", "--guess").help("-- Initial wavefunction.").default_value(std::string("exp(-(x-1)^2)"));
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(1000).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-n", "--nstate").help("-- Number of states to calculate.").default_value(1).scan<'i', int>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(0.1).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    int iters = program.get<int>("-i"), nstate = program.get<int>("-n"); double mass = program.get<double>("-m"), step = program.get<double>("-s");

    // load the potential function and initialize wfn
    EigenMatrix<> U = Eigen::LoadMatrix("U.mat"); Wavefunction wfn(program.get<std::string>("-g"), U.leftCols(U.cols() - 1), mass);

    // normalize the wfn, extract the potential values from the last column
    wfn = wfn / wfn.norm(); U.col(0) = U.rightCols(1); U.conservativeResize(U.rows(), 1);

    // define the imaginary time kinetic and potential operators
    EigenMatrix<std::complex<double>> R = (-0.5 * U.array() * step).exp(), K = (-0.5 * wfn.getk().array().pow(2) * step / mass).exp();

    // define the energy and wfn vector
    EigenVector<> eps(nstate); std::vector<Wavefunction> states(nstate, wfn);

    // perform the dynamics for every state
    for (int i = 0; i < nstate; i++) {

        // assign the initial energy and wfn
        wfn = states.at(i); double E = wfn.energy(U);

        // propagate the current state
        for (int j = 0; j < iters; j++) {

            // copy the wavefunction and energy
            Wavefunction wfnp = wfn; double Ep = E;

            // propagate the wavefunction
            wfn = wfn.propagate(R, K);

            // subrtract the lower eigenstates
            for (int k = 0; k < i; k++) wfn = wfn - states.at(k) * wfn.overlap(states.at(k));

            // normalize the wavefunction and calculate the energy
            wfn = wfn / wfn.norm(); E = wfn.energy(U);

            // print the iteration info
            std::printf("%6d %20.14f %.2e %.2e %s\n", j + 1, E, std::abs(E - Ep), (wfn - wfnp).norm(), "");
        }

        // append to states and energy and print the new line
        states.at(i) = wfn; eps(i) = E; std::cout << std::endl;
    }

    // print the final energies
    std::cout << eps << std::endl;
}

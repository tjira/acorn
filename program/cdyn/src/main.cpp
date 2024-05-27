#include "expression.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Classical Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-d", "--derivative").help("-- The the derivatives of the potentials in diabatic representation.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<std::string>{"x"});
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-c", "--coordinate").help("-- Initial position of the system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-r", "--random").help("-- Random seed for the simulation.").default_value(1).scan<'i', int>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("-t", "--trajectories").help("-- Number of trajectories to run.").default_value(1).scan<'i', int>();
    program.add_argument("-u", "--potential").help("-- The potential curves in diabatic representation.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<std::string>{"0.5*x^2"});
    program.add_argument("--savetraj").help("-- Save the trajectories.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters and define the potential expressions
    double momentum = program.get<double>("-p"), position = program.get<double>("-c"), step = program.get<double>("-s"), mass = program.get<double>("-m"); std::vector<Expression> U, dU;
    int iters = program.get<int>("-i"), seed = program.get<int>("-r"), trajs = program.get<int>("-t"); int nstate = std::sqrt(program.get<std::vector<std::string>>("-u").size());

    // extract boolean flags
    bool savetraj = program.get<bool>("--savetraj");

    // reserve the memory for the potential expressions
    U.reserve(nstate * nstate), dU.reserve(nstate * nstate);

    // fill the potential and potential derivative expressions
    for (const std::string& expr : program.get<std::vector<std::string>>("-d")) dU.emplace_back(expr, std::vector<std::string>{"x"});
    for (const std::string& expr : program.get<std::vector<std::string>>("-u")) U.emplace_back(expr, std::vector<std::string>{"x"});

    // function to obtain the potential matrices
    auto eval = [&](std::vector<Expression>& Ue, double r) {
        EigenMatrix<> U(nstate, nstate); for (int i = 0; i < nstate; i++) {for (int j = 0; j < nstate; j++) U(i, j) = Ue.at(i * nstate + j).eval(r);} return U;
    };

    // define the container for all the trajectories and state
    EigenMatrix<> rt, vt, at; EigenMatrix<int> st(iters + 1, trajs);

    // initialize the container for all the trajectories
    if (savetraj) rt = EigenMatrix<>(iters + 1, trajs), vt = EigenMatrix<>(iters + 1, trajs), at = EigenMatrix<>(iters + 1, trajs);

    // print the header
    std::printf("%6s %6s %6s %14s %14s %14s %14s %s\n", "TRAJ", "ITER", "STATE", "EPOT", "EKIN", "ETOT", "FORCE", "");

    // loop over all trajectories
    for (int i = 0; i < trajs; i++) {

        // initialize random number generators
        std::mt19937 mt(seed + i); std::uniform_real_distribution<double> dist(0, 1);
    
        // initialize the containers for the position, velocity, acceleration, state energy difference and state
        EigenMatrix<> r(iters + 1, 1), v(iters + 1, 1), a(iters + 1, 1), de(iters + 1, 1); EigenMatrix<int> s(iters + 1, 1);

        // distributions for position and momentum
        std::normal_distribution<double> positiondist(position, 0.5);
        std::normal_distribution<double> momentumdist(momentum, 1.0);

        // fill the initial position, velocity, acceleration state and state difference
        r(0) = positiondist(mt), v(0) = momentumdist(mt) / mass, a(0) = 0, s(0) = 0; if (nstate > 1) de(0) = eval(U, r(0))(1, 1) - eval(U, r(0))(0, 0);

        // calculate the force and potential with kinetic energy
        double F = -eval(dU, r(0))(s(0), s(0)), Epot = eval(U, r(0))(s(0), s(0)), Ekin = 0.5 * mass * v(0) * v(0);

        // print the zeroth iteration
        std::printf("%6d %6d %6d %14.8f %14.8f %14.8f %14.8f %s\n", i + 1, 0, s(0), Epot, Ekin, Epot + Ekin, F, "");

        // perform the dynamics
        for (int j = 0; j < iters; j++) {

            // start the timer
            if (!((j + 1) % 1)) tp = Timer::Now();

            // calculate the velocity and accceleration
            a(j + 1) = F / mass; v(j + 1) = v(j) + 0.5 * (a(j) + a(j + 1)) * step;

            // move the system and propagate state
            r(j + 1) = r(j) + step * (v(j + 1) + 0.5 * a(j + 1) * step), s(j + 1) = s(j);

            // nonadiabatic dynamics
            if (nstate > 1) {

                // calculate and save the difference in state energies
                de(j + 1) = eval(U, r(j + 1))(1, 1) - eval(U, r(j + 1))(0, 0);

                // calculate the probability of state change according to the Landau-Zener formula
                double gamma = std::pow(eval(U, r(j + 1))(0, 1), 2) / std::abs((de(j + 1) - de(j)) / step); double P = 1 - std::exp(-2 * M_PI * gamma);

                // change the state if the jump is accepted
                if (de(j) * de(j + 1) < 0 && dist(mt) < P) {
                    s(j + 1) = s(j + 1) == 1 ? 0 : 1; v(0) = std::sqrt(v(0)*v(0) - (s(j + 1) - s(j)) * 2 * de(j + 1) / mass);
                }
            }

            // calculate the potential and kinetic energy
            Epot = eval(U, r(j + 1))(s(j + 1), s(j + 1)), Ekin = 0.5 * mass * v(j + 1) * v(j + 1);

            // calculate the force
            F = -eval(dU, r(j + 1))(s(j + 1), s(j + 1));

            // print the iteration
            if (!((j + 1) % 1)) std::printf("%6d %6d %6d %14.8f %14.8f %14.8f %14.8f %s\n", i + 1, j + 1, s(j), Epot, Ekin, Epot + Ekin, F, "");
        }

        // save the trajectory to the containers
        if (savetraj) {rt.col(i) = r, vt.col(i) = v, at.col(i) = a;} st.col(i) = s;
    }

    // print the new line
    std::cout << std::endl;

    // create the container for the state populations
    EigenMatrix<> P(iters + 1, nstate + 1); P.setZero();

    // fill the state populations
    for (int i = 0; i < iters + 1; i++) {
        for (int j = 0; j < trajs; j++) {
            P(i, 1 + st(i, j)) += 1.0 / trajs;
        } P(i, 0) = i * step;
    }

    // write the trajectories to the disk
    MEASURE("TRAJECTORY AND POPULATION WRITING: ", 
        if (savetraj) Eigen::Write("TRAJ_POS.mat", rt);
        if (savetraj) Eigen::Write("TRAJ_VEL.mat", vt);
        if (savetraj) Eigen::Write("TRAJ_ACC.mat", at);
                     Eigen::Write("P.mat"       , P );
    )

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

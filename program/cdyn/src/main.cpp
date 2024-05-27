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
    program.add_argument("--printint").help("-- Printing interval of the simulations.").default_value(1).scan<'i', int>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters and define the potential expressions
    double momentum = program.get<double>("-p"), position = program.get<double>("-c"), step = program.get<double>("-s"), mass = program.get<double>("-m");
    int iters = program.get<int>("-i"), pi = program.get<int>("--printint"), seed = program.get<int>("-r"), trajs = program.get<int>("-t");
    int nstate = std::sqrt(program.get<std::vector<std::string>>("-u").size()); std::vector<Expression> U, dU;

    // fill the potential and potential derivative expressions
    for (const std::string& expr : program.get<std::vector<std::string>>("-d")) dU.push_back(Expression(expr, {"x"}));
    for (const std::string& expr : program.get<std::vector<std::string>>("-u")) U.push_back(Expression(expr, {"x"}));

    // function to obtain the potential matrices
    auto eval = [&](std::vector<Expression>& Ue, double r) {
        EigenMatrix<> U(nstate, nstate); for (int i = 0; i < nstate; i++) {for (int j = 0; j < nstate; j++) U(i, j) = Ue.at(i * nstate + j).eval(r);} return U;
    };

    // define the container for all the trajectories
    EigenMatrix<> rt(iters + 1, trajs), vt(iters + 1, trajs), at(iters + 1, trajs);

    // print the header
    std::printf("%6s %6s %14s %14s %14s %14s %s\n", "TRAJ", "ITER", "EPOT", "EKIN", "ETOT", "FORCE", "");

    // loop over all trajectories
    for (int i = 0; i < trajs; i++) {

        // random number generators
        std::mt19937 mt(seed + i); std::uniform_real_distribution<double> dist(0, 1);
    
        // create the containers for the position, velocity and acceleration
        EigenMatrix<> r(iters + 1, 1), v(iters + 1, 1), a(iters + 1, 1);

        // distributions for position and momentum
        std::normal_distribution<double> positiondist(position, 0.5);
        std::normal_distribution<double> momentumdist(momentum, 1.0);

        // fill the initial position and velocity
        r(0) = positiondist(mt), v(0) = momentumdist(mt) / mass;

        // calculate the force and potential with kinetic energy
        double F = -eval(dU, r(0))(0, 0), Epot = eval(U, r(0))(0, 0), Ekin = 0.5 * mass * v(0) * v(0);

        // print the zeroth iteration
        std::printf("%6d %6d %14.8f %14.8f %14.8f %14.8f %s\n", i + 1, 0, Epot, Ekin, Epot + Ekin, F, "");

        // perform the dynamics
        for (int j = 0; j < iters; j++) {

            // start the timer
            if (!((j + 1) % pi)) tp = Timer::Now();

            // calculate the velocity and accceleration
            a(j + 1) = F / mass; v(j + 1) = v(j) + 0.5 * (a(j) + a(j + 1)) * step;

            // move the system
            r(j + 1) = r(j) + step * (v(j + 1) + 0.5 * a(j + 1) * step);

            // calculate the potential and kinetic energy and define the hopping probability variable
            Epot = eval(U, r(j + 1))(0, 0), Ekin = 0.5 * mass * v(j + 1) * v(j + 1);

            // calculate the force
            F = -eval(dU, r(j + 1))(0, 0);

            // print the iteration
            if (!((j + 1) % pi)) std::printf("%6d %6d %14.8f %14.8f %14.8f %14.8f %s\n", i + 1, j + 1, Epot, Ekin, Epot + Ekin, F, "");
        }

        // save the trajectory to the containers
        rt.col(i) = r, vt.col(i) = v, at.col(i) = a;
    }

    // print the new line
    std::cout << std::endl;

    // write the trajectories to the disk
    MEASURE("TRAJECTORY WRITING: ", 
        Eigen::Write("TRAJ_POS.mat", rt);
        Eigen::Write("TRAJ_VEL.mat", vt);
        Eigen::Write("TRAJ_ACC.mat", at);
    )

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

#include "expression.h"
#include "timer.h"
#include <argparse.hpp>

std::vector<std::vector<int>> Combinations(int n, int k) {
    // create the bitmask that will get permuted and the resulting vector
    std::string bitmask(k, 1); bitmask.resize(n, 0); std::vector<std::vector<int>> combs;
 
    // generate the combinations
    do {std::vector<int> comb; comb.reserve(k);
        for (int j = 0; j < n; j++) {
            if (bitmask[j]) comb.push_back(j);
        } combs.push_back(comb);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    // return the result
    return combs;
}

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Classical Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-c", "--coordinate").help("-- Initial position of the system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-e", "--excstate").help("-- Initial state of the system.").default_value(0).scan<'i', int>();
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-r", "--random").help("-- Random seed for the simulation.").default_value(1).scan<'i', int>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(1.0).scan<'g', double>();
    program.add_argument("-t", "--trajectories").help("-- Number of trajectories to run.").default_value(1).scan<'i', int>();
    program.add_argument("-u", "--potential").help("-- The potential curves in diabatic representation.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<std::string>{"0.5*x^2"});
    program.add_argument("--adiabatic").help("-- Perform the dynamics in adiabatic basis.").default_value(false).implicit_value(true);
    program.add_argument("--savetraj").help("-- Save the trajectories.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // extract the command line parameters and define the potential expressions
    double momentum = program.get<double>("-p"), position = program.get<double>("-c"), step = program.get<double>("-s"), mass = program.get<double>("-m"); std::vector<Expression> Ue, dUe;
    int iters = program.get<int>("-i"), seed = program.get<int>("-r"), trajs = program.get<int>("-t"); int nstate = std::sqrt(program.get<std::vector<std::string>>("-u").size());

    // extract boolean flags
    bool adiabatic = program.get<bool>("--adiabatic"), savetraj = program.get<bool>("--savetraj");

    // reserve the memory for the potential expressions
    Ue.reserve(nstate * nstate), dUe.reserve(nstate * nstate);

    // fill the potential expressions
    for (const std::string& expr : program.get<std::vector<std::string>>("-u")) Ue.emplace_back(expr, std::vector<std::string>{"x"});

    // function to obtain the potential matrices
    auto eval = [&](std::vector<Expression>& Ue, double r) {
        Matrix U(nstate, nstate); for (int i = 0; i < nstate; i++) {for (int j = 0; j < nstate; j++) U(i, j) = Ue.at(i * nstate + j).eval(r);} return U;
    };

    // define the container for all the trajectories and state
    Matrix rt, vt, at; IntegerMatrix st(iters + 1, trajs);

    // initialize the container for all the trajectories
    if (savetraj) rt = Matrix(iters + 1, trajs + 1);

    // print the header
    std::printf("%6s %6s %6s %14s %14s %14s %14s %s\n", "TRAJ", "ITER", "STATE", "EPOT", "EKIN", "ETOT", "FORCE", "");

    // loop over all trajectories
    for (int i = 0; i < trajs; i++) {

        // initialize random number generators and the potential matrix with its derivative
        std::mt19937 mt(seed + i); std::uniform_real_distribution<double> dist(0, 1); Matrix U, dU;
    
        // initialize the containers for the position, velocity, acceleration, state energy difference and state
        Matrix r(iters + 1, 1), v(iters + 1, 1), a(iters + 1, 1), de(iters + 1, Combinations(nstate, 2).size()), dz(iters + 1, Combinations(nstate, 2).size()); IntegerMatrix s(iters + 1, 1);

        // distributions for position and momentum
        std::normal_distribution<double> positiondist(position, 0.5);
        std::normal_distribution<double> momentumdist(momentum, 1.0);

        // fill the initial position, velocity, acceleration and state
        r(0) = positiondist(mt), v(0) = momentumdist(mt) / mass, a(0) = 0, s(0) = program.get<int>("-e");

        // calculate the potential and perform adiabatic transform
        if (U = eval(Ue, r(0)); adiabatic) {Eigen::SelfAdjointEigenSolver<Matrix> solver(U); U = solver.eigenvalues().asDiagonal();};

        // define the containers for the potential, save the first value and define the state combinations
        std::vector<Matrix> Uc(iters + 1); Uc.at(0) = U; auto scombs = Combinations(nstate, 2); 

        // calcualte the state energy difference
        if (nstate > 1) for (int j = 0; j < de.cols(); j++) {
            de(0, j) = U(scombs.at(j).at(1), scombs.at(j).at(1)) - U(scombs.at(j).at(0), scombs.at(j).at(0));
        }

        // initialize the force and potential with kinetic energy
        double F = mass * a(0), Epot = U(s(0), s(0)), Ekin = 0.5 * mass * v(0) * v(0);

        // save the initial point in phase space
        if (savetraj) {rt(0, 0) = 0, rt(0, i + 1) = Epot;}

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

            // calculate the potential and perform adiabatic transform
            if (U = eval(Ue, r(j + 1)); adiabatic) {Eigen::SelfAdjointEigenSolver<Matrix> solver(U); U = solver.eigenvalues().asDiagonal();};

            // calculate the gradient of the potential and save U and dU to the container
            Uc.at(j + 1) = U, dU = (Uc.at(j + 1) - Uc.at(j)) / (r(j + 1) - r(j));

            // calculate state energy difference and its numerical derivative
            if (nstate > 1) for (int k = 0; k < de.cols(); k++) {
                de(j + 1, k) = U(scombs.at(k).at(1), scombs.at(k).at(1)) - U(scombs.at(k).at(0), scombs.at(k).at(0)); dz(j + 1, k) = (de(j + 1, k) - de(j, k)) / step;
            }

            // nonadiabatic dynamics
            if (double P; nstate > 1 && j > 0) {

                // for every state cobination
                for (int k = 0; k < de.cols(); k++) {

                    // skip if the current state is not in the combination
                    if (s(j + 1) != scombs.at(k).at(0) && s(j + 1) != scombs.at(k).at(1)) continue;

                    // calculate the probability of state change according to the Landau-Zener formula
                    if (adiabatic) {
                        double ddz = (de(j + 1, k) - 2 * de(j, k) + de(j - 1, k)) / step / step; P = std::exp(-0.5 * M_PI * std::sqrt(std::pow(de(j + 1, k), 3) / ddz));
                    } else {
                        double gamma = std::pow(U(scombs.at(k).at(0), scombs.at(k).at(1)), 2) / std::abs((de(j + 1, k) - de(j, k)) / step); P = 1 - std::exp(-2 * M_PI * gamma);
                    }

                    // change the state and rescale the velocities if the jump is accepted
                    if ((adiabatic && dz(j, k) * dz(j + 1, k) < 0 && dist(mt) < P) || (!adiabatic && de(j, k) * de(j + 1, k) < 0 && dist(mt) < P)) {

                        // change the state
                        s(j + 1) = s(j + 1) == scombs.at(k).at(0) ? scombs.at(k).at(1) : scombs.at(k).at(0);

                        // rescale the velocity
                        int factor = s(j + 1) > s(j) ? 1 : -1; v(0) = std::sqrt(v(0) * v(0) - factor * 2 * de(j + 1, k) / mass);
                    }
                }
            }

            // calculate the potential and kinetic energy
            Epot = U(s(j + 1), s(j + 1)), Ekin = 0.5 * mass * v(j + 1) * v(j + 1);

            // calculate the force
            F = -dU(s(j + 1), s(j + 1));

            // save the point of trajectory
            if (savetraj) rt(j + 1, 0) = (j + 1) * step, rt(j + 1, i + 1) = Epot;

            // print the iteration
            if (!((j + 1) % 1)) std::printf("%6d %6d %6d %14.8f %14.8f %14.8f %14.8f %s\n", i + 1, j + 1, s(j), Epot, Ekin, Epot + Ekin, F, "");
        }

        // assign state
        st.col(i) = s;
    }

    // create the container for the state populations
    Matrix P(iters + 1, nstate + 1); P.setZero();

    // fill the state populations
    for (int i = 0; i < iters + 1; i++) {
        for (int j = 0; j < trajs; j++) {
            P(i, 1 + st(i, j)) += 1.0 / trajs;
        } P(i, 0) = i * step;
    }

    // print the new line
    std::cout << std::endl;

    // print the populations
    for (int i = 0; i < P.cols() - 1; i++) {
        std::printf("%s STATE %d POP: %.14f\n", adiabatic ? "ADIABATIC" : "DIABATIC", i, P(P.rows() - 1, i + 1));
    } std::cout << std::endl;

    // write the trajectories to the disk
    MEASURE("TRAJECTORY AND POPULATION WRITING: ", 
        if (savetraj) {Eigen::Write(std::string("TRAJ_LZ_") + (adiabatic ? "ADIA" : "DIA") + ".mat", rt);} Eigen::Write(std::string("P_LZ_") + (adiabatic ? "ADIA" : "DIA") + ".mat", P);
    )

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

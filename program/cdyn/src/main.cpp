#include "expression.h"
#include "landauzener.h"
#include "timer.h"
#include <argparse.hpp>
#include <format>

int nthread = 1;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Classical Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-c", "--coordinate").help("-- Initial position of the system.").default_value(0.0).scan<'g', double>();
    program.add_argument("-e", "--excstate").help("-- Initial diabatic state of the system.").default_value(0).scan<'i', int>();
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(200).scan<'i', int>();
    program.add_argument("-l", "--log").help("-- Log interval for the simulation.").default_value(10).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(1.0).scan<'g', double>();
    program.add_argument("-n", "--nthreads").help("-- Number of threads to use.").default_value(nthread).scan<'i', int>();
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

    // extract the command line parameters
    double momentum = program.get<double>("-p"), position = program.get<double>("-c"), step = program.get<double>("-s"), mass = program.get<double>("-m"); nthread = program.get<int>("-n");
    int iters = program.get<int>("-i"), seed = program.get<int>("-r"), trajs = program.get<int>("-t"); int nstate = std::sqrt(program.get<std::vector<std::string>>("-u").size());

    // extract boolean flags and the log interval
    bool adiabatic = program.get<bool>("--adiabatic"), savetraj = program.get<bool>("--savetraj"); int log = program.get<int>("-l");

    // define the potential expressions and allocate the memory
    std::vector<std::vector<Expression>> Ues(nthread); for (auto& Ue : Ues) Ue.reserve(nstate * nstate);

    // fill the potential expressions
    for (auto& Ue : Ues) for (const std::string& expr : program.get<std::vector<std::string>>("-u")) Ue.emplace_back(expr, std::vector<std::string>{"x"});

    // function to obtain the potential matrices
    auto eval = [&](std::vector<Expression>& Ue, double r) {
        Matrix U(nstate, nstate); for (int i = 0; i < nstate; i++) {for (int j = 0; j < nstate; j++) U(i, j) = Ue.at(i * nstate + j).eval(r);} return U;
    };

    // define the container for all the trajectories and states
    Matrix sdt(iters + 1, nstate * trajs), sat; if (adiabatic) sat = sdt;

    // print the header
    std::printf("%6s %6s %6s %14s %14s %14s %14s\n", "TRAJ", "ITER", "STATE", "EPOT", "EKIN", "ETOT", "FORCE");

    // loop over all trajectories
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < trajs; i++) {

        // get the thread id
        #ifdef _OPENMP
        int id = omp_get_thread_num();
        #else
        int id = 0;
        #endif

        // distributions for position and momentum, mersenne twister and uniform distribution
        std::mt19937 mt(trajs * (seed + i) + 1); std::uniform_real_distribution<double> dist(0, 1);
        std::normal_distribution<double> positiondist(position, 0.5), momentumdist(momentum, 1.0);

        // initialize the containers for the position, velocity, acceleration, state energy difference, state and a random number
        Matrix U, dU, r(iters + 1, 1), v(iters + 1, 1), a(iters + 1, 1); IntegerMatrix sd(iters + 1, 1); double rn;

        // fill the initial position, velocity, acceleration and state and evaluate the potential
        r(0) = positiondist(mt), v(0) = momentumdist(mt) / mass, a(0) = 0, sd(0) = program.get<int>("-e"); U = eval(Ues.at(id), r(0));

        // initialize the Landau-Zener algorithm and the transition probabilities
        LandauZener lz = nstate > 1 ? LandauZener(nstate, iters + 1) : LandauZener(); std::vector<std::tuple<int, double, bool>> transitions;

        // define the container for the potential and attempt a LZ jump
        std::vector<Matrix> Uc(iters + 1); Uc.at(0) = U; if (nstate > 1) transitions = lz.jump(U, sd(0), 0, step);

        // initialize the force and potential with kinetic energy
        double F = mass * a(0), Epot = U(sd(0), sd(0)), Ekin = 0.5 * mass * v(0) * v(0);

        // print the zeroth iteration
        std::printf("%6d %6d %6d %14.8f %14.8f %14.8f %14.8f\n", i + 1, 0, sd(0), Epot, Ekin, Epot + Ekin, F);

        // perform the dynamics
        for (int j = 0; j < iters; j++) {

            // start the timer
            if (!((j + 1) % 1)) tp = Timer::Now();

            // calculate the velocity and accceleration
            a(j + 1) = F / mass; v(j + 1) = v(j) + 0.5 * (a(j) + a(j + 1)) * step;

            // move the system and propagate state
            r(j + 1) = r(j) + step * (v(j + 1) + 0.5 * a(j + 1) * step), sd(j + 1) = sd(j);

            // evaluate the potential
            U = eval(Ues.at(id), r(j + 1)); Uc.at(j + 1) = U, dU = (Uc.at(j + 1) - Uc.at(j)) / (r(j + 1) - r(j));

            // calculate the Landau-Zener jump probabilities
            if (nstate > 1) transitions = lz.jump(U, sd(j + 1), j + 1, step);

            // loop over all the transitions
            if (rn = dist(mt); nstate > 1) for (const auto& transition : transitions) {
                if (int ns = std::get<0>(transition); std::get<2>(transition) && rn < std::get<1>(transition)) {
                    v(j + 1) = std::sqrt(v(j + 1) * v(j + 1) - 2 * (U(ns, ns) - U(sd(j + 1), sd(j + 1))) / mass); sd(j + 1) = ns; break;
                }
            }

            // calculate the potential and kinetic energy
            Epot = U(sd(j + 1), sd(j + 1)), Ekin = 0.5 * mass * v(j + 1) * v(j + 1);

            // calculate the force
            F = -dU(sd(j + 1), sd(j + 1));

            // print the iteration
            if (std::string probs; !((j + 1) % log)) {
                for (size_t k = 0; k < transitions.size(); k++) {
                    probs += std::format("{}{}->{}={:.3f}({})", k ? ", " : " ", sd(j), std::get<0>(transitions.at(k)), std::get<1>(transitions.at(k)), std::get<2>(transitions.at(k)));
                }
                std::printf("%6d %6d %6d %14.8f %14.8f %14.8f %14.8f %.4f %s\n", i + 1, j + 1, sd(j), Epot, Ekin, Epot + Ekin, F, rn, probs.c_str());
            }
        }

        // assign the diabatic states to the time dependent container
        for (int j = 0; j < iters + 1; j++) for (int k = 0; k < nstate; k++) sdt(j, i * nstate + k) = k == sd(j) ? 1 : 0;

        // assign the adiabatic states to the time dependent container
        for (int j = 0; j < iters + 1 && adiabatic; j++) {

            // diagonalize the potential for each time step of the trajectory
            Eigen::SelfAdjointEigenSolver<Matrix> solver(Uc.at(j)); Matrix C = solver.eigenvectors();

            // assign the adiabatic state to the time dependent container
            sat.block(j, i * nstate, 1, nstate) = (C.adjoint() * sdt.block(j, i * nstate, 1, nstate).transpose().asDiagonal().toDenseMatrix() * C).diagonal().transpose();
        }

        // save trajectories and all energy data
        if (savetraj) {

            // create the matrices
            Matrix E(iters + 1, 2), ED(iters + 1, lz.getEd().cols() + 1), DED(iters + 1, lz.getDed().cols() + 1), RT(iters + 1, 2), VT(iters + 1, 2), AT(iters + 1, 2);

            // fill the position, velocity and acceleration matrices
            RT << Vector::LinSpaced(iters + 1, 0, iters) * step, r; VT << Vector::LinSpaced(iters + 1, 0, iters) * step, v; AT << Vector::LinSpaced(iters + 1, 0, iters) * step, a;

            // write the position, velocity and acceleration matrices to disk
            Eigen::Write("R_LZ_" + std::to_string(i + 1) + ".mat", RT); Eigen::Write("V_LZ_" + std::to_string(i + 1) + ".mat", VT); Eigen::Write("A_LZ_" + std::to_string(i + 1) + ".mat", AT);

            // fill the energy matrix
            for (int j = 0; j < iters + 1; j++) E(j, 0) = j * step, E(j, 1) = Uc.at(j)(sd(j), sd(j));

            // fill the energy difference matrices
            ED << Vector::LinSpaced(iters + 1, 0, iters) * step, lz.getEd(); DED << Vector::LinSpaced(iters + 1, 0, iters) * step, lz.getDed();

            // write the energy matrices to disk
            Eigen::Write("E_LZ_" + std::to_string(i + 1) + ".mat", E); Eigen::Write("ED_LZ_" + std::to_string(i + 1) + ".mat", ED); Eigen::Write("DED_LZ_" + std::to_string(i + 1) + ".mat", DED);
        }
    }

    // create the container for the diabatic and adiabatic state populations and fill the time
    Matrix Pd(iters + 1, nstate * nstate + 1), Pa; Pd.setZero(); for (int i = 0; i < iters + 1; i++) Pd(i, 0) = i * step; if (adiabatic) Pa = Pd;

    // fill the diabatic state populations
    for (int i = 0; i < nstate; i++) {
        Pd.col(1 + i * (nstate + 1)) = sdt(Eigen::all, Vector::LinSpaced(trajs, 0, trajs - 1) * nstate + Vector::Constant(trajs, i)).rowwise().sum() / trajs;
    }

    // fill the adiabatic state populations
    for (int i = 0; i < nstate && adiabatic; i++) {
        Pa.col(1 + i * (nstate + 1)) = sat(Eigen::all, Vector::LinSpaced(trajs, 0, trajs - 1) * nstate + Vector::Constant(trajs, i)).rowwise().sum() / trajs;
    }

    // print the new line
    std::cout << std::endl;

    // print the diabatic populations
    for (int i = 0; i < nstate; i++) {
        std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(nstate * nstate).reshaped(nstate, nstate)(i, i));
    } std::cout << std::endl;

    // print the adiabatic populations
    for (int i = 0; i < nstate && adiabatic; i++) {
        std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(nstate * nstate).reshaped(nstate, nstate)(i, i));
    } std::cout << std::endl;

    // write the trajectories to disk
    MEASURE("POPULATION WRITING: ", 
        Eigen::Write("P_LZ_DIA.mat", Pd); if (adiabatic) Eigen::Write("P_LZ_ADIA.mat", Pa);
    )

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

#include "cdyn.h"
#include "expression.h"
#include "lz.h"
#include "timer.h"
#include <argparse.hpp>

int nthread = 1;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Classical Dynamics Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

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
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set the number of threads end extract the number of states
    nthread = program.get<int>("-n"); int nstate = std::sqrt(program.get<std::vector<std::string>>("-u").size());

    // extract the command line parameters
    Acorn::CDYN::Options opt = {
        program.get<std::vector<std::string>>("-u"), program.get<double>("-m"         ), program.get<double>("-p"        ), program.get<double>("-c"), program.get<double>("-s"),
        program.get<int                     >("-e"), program.get<int   >("-i"         ), program.get<int   >("-l"        ), nstate,                    program.get<int   >("-r"),
        program.get<int                     >("-t"), program.get<bool  >("--adiabatic"), program.get<bool  >("--savetraj")
    };

    // define the potential expressions and allocate the memory
    std::vector<std::vector<Expression>> Ues(nthread); for (auto& Ue : Ues) Ue.reserve(nstate * nstate);

    // fill the potential expressions
    for (auto& Ue : Ues) for (const std::string& expr : opt.potential) Ue.emplace_back(expr, std::vector<std::string>{"x"});

    // function to obtain the potential matrices
    auto eval = [&](std::vector<Expression>& Ue, double r) {
        Matrix U(nstate, nstate); for (int i = 0; i < nstate; i++) {for (int j = 0; j < nstate; j++) U(i, j) = Ue.at(i * nstate + j).eval(r);} return U;
    };

    // define the container for all the states and hopping geometries
    Matrix sdt(opt.iters + 1, nstate * opt.trajs), sat(opt.iters + 1, nstate * opt.trajs); std::vector<double> hg;

    // print the header
    std::printf("%6s %6s %6s %14s %14s %14s %14s\n", "TRAJ", "ITER", "STATE", "EPOT", "EKIN", "ETOT", "FORCE");

    // loop over all trajectories
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < opt.trajs; i++) {

        // get the thread id
        #ifdef _OPENMP
        int id = omp_get_thread_num();
        #else
        int id = 0;
        #endif

        // distributions for position and momentum, mersenne twister and uniform distribution
        std::mt19937 mt(opt.trajs * (opt.seed + i) + 1); std::uniform_real_distribution<double> dist(0, 1);
        std::normal_distribution<double> positiondist(opt.position, 0.5), momentumdist(opt.momentum, 1.0);

        // initialize the containers for the position, velocity, acceleration, state energy difference, state and a random number
        Matrix U, dU, r(opt.iters + 1, 1), v(opt.iters + 1, 1), a(opt.iters + 1, 1); IntegerMatrix s(opt.iters + 1, 1); double rn;

        // fill the initial position, velocity, acceleration and state
        r(0) = positiondist(mt), v(0) = momentumdist(mt) / opt.mass, a(0) = 0, s(0) = program.get<int>("-e");

        // diagonalize the potential at the initial position
        U = eval(Ues.at(id), r(0)); if (opt.adiabatic) {Eigen::SelfAdjointEigenSolver<Matrix> solver(U); U = solver.eigenvalues().asDiagonal();}

        // initialize the Landau-Zener algorithm and the transition probabilities
        LandauZener lz = nstate > 1 ? LandauZener(nstate, opt.iters + 1, opt.adiabatic) : LandauZener(); std::vector<std::tuple<int, double, bool>> transitions;

        // define the container for the potential and attempt a LZ jump
        std::vector<Matrix> Uc(opt.iters + 1); Uc.at(0) = U; if (nstate > 1) transitions = lz.jump(U, s(0), 0, opt.step);

        // initialize the force and potential with kinetic energy
        double F = opt.mass * a(0), Epot = U(s(0), s(0)), Ekin = 0.5 * opt.mass * v(0) * v(0);

        // print the zeroth iteration
        std::printf("%6d %6d %6d %14.8f %14.8f %14.8f %14.8f\n", i + 1, 0, s(0), Epot, Ekin, Epot + Ekin, F);

        // perform the dynamics
        for (int j = 0; j < opt.iters; j++) {

            // calculate the velocity and accceleration
            a(j + 1) = F / opt.mass; v(j + 1) = v(j) + 0.5 * (a(j) + a(j + 1)) * opt.step;

            // move the system and set the next state to the same as before
            r(j + 1) = r(j) + opt.step * (v(j + 1) + 0.5 * a(j + 1) * opt.step), s(j + 1) = s(j);

            // evaluate the potential and transform it to adiabatic basis if needed
            U = eval(Ues.at(id), r(j + 1)); if (opt.adiabatic) {Eigen::SelfAdjointEigenSolver<Matrix> solver(U); U = solver.eigenvalues().asDiagonal();}

            // store the potential and calculate its derivative
            Uc.at(j + 1) = U, dU = (Uc.at(j + 1) - Uc.at(j)) / (r(j + 1) - r(j));

            // calculate the Landau-Zener jump probabilities
            if (nstate > 1) transitions = lz.jump(U, s(j + 1), j + 1, opt.step);

            // loop over all the transitions
            if (rn = dist(mt); nstate > 1) for (size_t k = 0; k < transitions.size(); k++) { // loop over all the transitions
                if (int ns = std::get<0>(transitions.at(k)); std::get<2>(transitions.at(k))) { // if the jump attempt is performed
                    if (rn > (k ? std::get<1>(transitions.at(k - 1)) : 0) && rn < std::get<1>(transitions.at(k))) { // if the jump is accepted
                        if (double vn = v(j + 1) * v(j + 1) - 2 * (U(ns, ns) - U(s(j + 1), s(j + 1))) / opt.mass; vn > 0) { // if the system has enough kinetic energy
                            v(j + 1) = std::sqrt(v(j + 1) * v(j + 1) - 2 * (U(ns, ns) - U(s(j + 1), s(j + 1))) / opt.mass); s(j + 1) = ns; hg.push_back(r(j + 1)); break;
                        }
                    }
                }
            }

            // calculate the potential and kinetic energy
            Epot = U(s(j + 1), s(j + 1)), Ekin = 0.5 * opt.mass * v(j + 1) * v(j + 1);

            // calculate the force
            F = -dU(s(j + 1), s(j + 1));

            // print the iteration
            if (std::string print; !((j + 1) % opt.log) || s(j+1) != s(j)) {
                char buffer[100]; std::sprintf(buffer, "%6d %6d %6d %14.8f %14.8f %14.8f %14.8f %.4f", i + 1, j + 1, s(j), Epot, Ekin, Epot + Ekin, F, rn); print = buffer;
                for (size_t k = 0; k < transitions.size(); k++) {
                    std::sprintf(buffer, "%s%d->%d=%.3f(%d)", k ? ", " : " ", s(j), std::get<0>(transitions.at(k)), std::get<1>(transitions.at(k)), std::get<2>(transitions.at(k))); print += buffer;
                } std::printf("%s\n", print.c_str());
            }
        }

        // assign the diabatic states to the time dependent container
        for (int j = 0; j < opt.iters + 1 && !opt.adiabatic; j++) for (int k = 0; k < nstate; k++) sdt(j, i * nstate + k) = k == s(j) ? 1 : 0;
        for (int j = 0; j < opt.iters + 1 &&  opt.adiabatic; j++) for (int k = 0; k < nstate; k++) sat(j, i * nstate + k) = k == s(j) ? 1 : 0;

        // assign the states to the time dependent container
        for (int j = 0; j < opt.iters + 1; j++) {

            // diagonalize the potential for each time step of the trajectory
            Eigen::SelfAdjointEigenSolver<Matrix> solver(eval(Ues.at(id), r(j))); Matrix C = solver.eigenvectors();

            // assign the state to the time dependent container
            if (!opt.adiabatic) sat.block(j, i * nstate, 1, nstate) = (C.adjoint() * sdt.block(j, i * nstate, 1, nstate).transpose().asDiagonal().toDenseMatrix() * C).diagonal().transpose();
            if ( opt.adiabatic) sdt.block(j, i * nstate, 1, nstate) = (C * sat.block(j, i * nstate, 1, nstate).transpose().asDiagonal().toDenseMatrix() * C.adjoint()).diagonal().transpose();
        }

        // save trajectories and all energy data
        if (opt.savetraj) {

            // create the energy
            Matrix E(opt.iters + 1, 2), ETOT(opt.iters + 1, 2), ED(opt.iters + 1, lz.getEd().cols() + 1), DED(opt.iters + 1, lz.getDed().cols() + 1), DDED(opt.iters + 1, lz.getDded().cols() + 1);

            // create the position, velocity and acceleration matrices
            Matrix RT(opt.iters + 1, 2), VT(opt.iters + 1, 2), AT(opt.iters + 1, 2);

            // fill the position, velocity and acceleration matrices
            RT << Vector::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, r;
            VT << Vector::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, v;
            AT << Vector::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, a;

            // write the position, velocity and acceleration matrices to disk
            Eigen::Write(std::string("R_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", RT);
            Eigen::Write(std::string("V_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", VT);
            Eigen::Write(std::string("A_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", AT);

            // fill the energy matrix
            for (int j = 0; j < opt.iters + 1; j++) E(j, 0) = j * opt.step, ETOT(j, 0) = j * opt.step, E(j, 1) = Uc.at(j)(s(j), s(j)), ETOT(j, 1) = Uc.at(j)(s(j), s(j)) + 0.5 * opt.mass * v(j) * v(j);

            // fill the energy difference matrices
            DDED << Vector::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, lz.getDded();
            DED  << Vector::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, lz.getDed();
            ED   << Vector::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, lz.getEd();

            // write the potential energy matrices to disk
            Eigen::Write(std::string("DDED_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", DDED);
            Eigen::Write(std::string("DED_LZ-" ) + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", DED );
            Eigen::Write(std::string("ED_LZ-"  ) + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", ED  );
            Eigen::Write(std::string("E_LZ-"   ) + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", E  );

            // write the total energy matrix to disk
            Eigen::Write(std::string("ETOT_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", ETOT);
        }
    }

    // create the container for the diabatic and adiabatic state populations and fill the time
    Matrix Pd(opt.iters + 1, nstate * nstate + 1), Pa; Pd.setZero(); for (int i = 0; i < opt.iters + 1; i++) Pd(i, 0) = i * opt.step; Pa = Pd;

    // fill the diabatic state populations
    for (int i = 0; i < nstate; i++) {
        Pd.col(1 + i * (nstate + 1)) = sdt(Eigen::all, Vector::LinSpaced(opt.trajs, 0, opt.trajs - 1) * nstate + Vector::Constant(opt.trajs, i)).rowwise().sum() / opt.trajs;
    }

    // fill the adiabatic state populations
    for (int i = 0; i < nstate; i++) {
        Pa.col(1 + i * (nstate + 1)) = sat(Eigen::all, Vector::LinSpaced(opt.trajs, 0, opt.trajs - 1) * nstate + Vector::Constant(opt.trajs, i)).rowwise().sum() / opt.trajs;
    }

    // print the new line
    std::cout << std::endl;

    // print the diabatic populations
    for (int i = 0; i < nstate; i++) {
        std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(nstate * nstate).reshaped(nstate, nstate)(i, i));
    } std::cout << std::endl;

    // print the adiabatic populations
    for (int i = 0; i < nstate; i++) {
        std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(nstate * nstate).reshaped(nstate, nstate)(i, i));
    } std::cout << std::endl;

    // write the populations to disk
    MEASURE("POPULATION WRITING: ", 
        Eigen::Write(std::string("P_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_DIA.mat", Pd); Eigen::Write(std::string("P_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_ADIA.mat", Pa);
    )

    // write the hopping geometries to disk
    Eigen::Write(std::string("HG_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_.mat", Eigen::Map<Matrix>(hg.data(), hg.size(), 1));

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}

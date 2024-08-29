#include "cdyn.h"

void Acorn::CDYN::run(const Options& opt, std::vector<timepoint>& timers) {
    // define the potential expressions and allocate the memory
    std::vector<std::vector<Expression>> Ues(Acorn::CDYN::nthread); for (auto& Ue : Ues) Ue.reserve(opt.nstate * opt.nstate);

    // fill the potential expressions
    for (auto& Ue : Ues) for (const std::string& expr : opt.potential) Ue.emplace_back(expr, std::vector<std::string>{"x"});

    // function to obtain the potential matrices
    auto eval = [&](std::vector<Expression>& Ue, double r) {
        Eigen::MatrixXd U(opt.nstate, opt.nstate); for (int i = 0; i < opt.nstate; i++) {for (int j = 0; j < opt.nstate; j++) U(i, j) = Ue.at(i * opt.nstate + j).eval(r);} return U;
    };

    // define the container for all the states and hopping geometries
    Eigen::MatrixXd sdt(opt.iters + 1, opt.nstate * opt.trajs), sat(opt.iters + 1, opt.nstate * opt.trajs); std::vector<double> hg;

    // print the header
    std::printf("%6s %6s %6s %14s %14s %14s %14s\n", "TRAJ", "ITER", "STATE", "EPOT", "EKIN", "ETOT", "FORCE");

    // loop over all trajectories
    #pragma omp parallel for num_threads(Acorn::CDYN::nthread)
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
        Eigen::MatrixXd U, dU, r(opt.iters + 1, 1), v(opt.iters + 1, 1), a(opt.iters + 1, 1); Eigen::MatrixXi s(opt.iters + 1, 1); double rn;

        // fill the initial position, velocity, acceleration and state
        r(0) = positiondist(mt), v(0) = momentumdist(mt) / opt.mass, a(0) = 0, s(0) = opt.excstate;

        // diagonalize the potential at the initial position
        U = eval(Ues.at(id), r(0)); if (opt.adiabatic) {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(U); U = solver.eigenvalues().asDiagonal();}

        // initialize the Landau-Zener algorithm and the transition probabilities
        Acorn::CDYN::LandauZener lz = opt.nstate > 1 ? LandauZener(opt.nstate, opt.iters + 1, opt.adiabatic) : Acorn::CDYN::LandauZener(); std::vector<std::tuple<int, double, bool>> transitions;

        // define the container for the potential and attempt a LZ jump
        std::vector<Eigen::MatrixXd> Uc(opt.iters + 1); Uc.at(0) = U; if (opt.nstate > 1) transitions = lz.jump(U, s(0), 0, opt.step);

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
            U = eval(Ues.at(id), r(j + 1)); if (opt.adiabatic) {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(U); U = solver.eigenvalues().asDiagonal();}

            // store the potential and calculate its derivative
            Uc.at(j + 1) = U, dU = (Uc.at(j + 1) - Uc.at(j)) / (r(j + 1) - r(j));

            // calculate the Landau-Zener jump probabilities
            if (opt.nstate > 1) transitions = lz.jump(U, s(j + 1), j + 1, opt.step);

            // loop over all the transitions
            if (rn = dist(mt); opt.nstate > 1) for (size_t k = 0; k < transitions.size(); k++) { // loop over all the transitions
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
        for (int j = 0; j < opt.iters + 1 && !opt.adiabatic; j++) for (int k = 0; k < opt.nstate; k++) sdt(j, i * opt.nstate + k) = k == s(j) ? 1 : 0;
        for (int j = 0; j < opt.iters + 1 &&  opt.adiabatic; j++) for (int k = 0; k < opt.nstate; k++) sat(j, i * opt.nstate + k) = k == s(j) ? 1 : 0;

        // assign the states to the time dependent container
        for (int j = 0; j < opt.iters + 1; j++) {

            // diagonalize the potential for each time step of the trajectory
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(eval(Ues.at(id), r(j))); Eigen::MatrixXd C = solver.eigenvectors();

            // assign the state to the time dependent container
            if (!opt.adiabatic) sat.block(j, i * opt.nstate, 1, opt.nstate) = (C.adjoint() * sdt.block(j, i * opt.nstate, 1, opt.nstate).transpose().asDiagonal().toDenseMatrix() * C).diagonal().transpose();
            if ( opt.adiabatic) sdt.block(j, i * opt.nstate, 1, opt.nstate) = (C * sat.block(j, i * opt.nstate, 1, opt.nstate).transpose().asDiagonal().toDenseMatrix() * C.adjoint()).diagonal().transpose();
        }

        // save trajectories and all energy data
        if (opt.savetraj) {

            // create the energy
            Eigen::MatrixXd E(opt.iters + 1, 2), ETOT(opt.iters + 1, 2), ED(opt.iters + 1, lz.getEd().cols() + 1), DED(opt.iters + 1, lz.getDed().cols() + 1), DDED(opt.iters + 1, lz.getDded().cols() + 1);

            // create the position, velocity and acceleration matrices
            Eigen::MatrixXd RT(opt.iters + 1, 2), VT(opt.iters + 1, 2), AT(opt.iters + 1, 2);

            // fill the position, velocity and acceleration matrices
            RT << Eigen::VectorXd::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, r;
            VT << Eigen::VectorXd::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, v;
            AT << Eigen::VectorXd::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, a;

            // write the position, velocity and acceleration matrices to disk
            torch::WriteMatrix(std::string("R_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", RT);
            torch::WriteMatrix(std::string("V_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", VT);
            torch::WriteMatrix(std::string("A_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", AT);

            // fill the energy matrix
            for (int j = 0; j < opt.iters + 1; j++) E(j, 0) = j * opt.step, ETOT(j, 0) = j * opt.step, E(j, 1) = Uc.at(j)(s(j), s(j)), ETOT(j, 1) = Uc.at(j)(s(j), s(j)) + 0.5 * opt.mass * v(j) * v(j);

            // fill the energy difference matrices
            DDED << Eigen::VectorXd::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, lz.getDded();
            DED  << Eigen::VectorXd::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, lz.getDed();
            ED   << Eigen::VectorXd::LinSpaced(opt.iters + 1, 0, opt.iters) * opt.step, lz.getEd();

            // write the potential energy matrices to disk
            torch::WriteMatrix(std::string("DDED_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", DDED);
            torch::WriteMatrix(std::string("DED_LZ-" ) + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", DED );
            torch::WriteMatrix(std::string("ED_LZ-"  ) + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", ED  );
            torch::WriteMatrix(std::string("E_LZ-"   ) + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", E  );

            // write the total energy matrix to disk
            torch::WriteMatrix(std::string("ETOT_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_" + std::to_string(i + 1) + ".mat", ETOT);
        }
    }

    // create the container for the diabatic and adiabatic state populations and fill the time
    Eigen::MatrixXd Pd(opt.iters + 1, opt.nstate * opt.nstate + 1), Pa; Pd.setZero(); for (int i = 0; i < opt.iters + 1; i++) Pd(i, 0) = i * opt.step; Pa = Pd;

    // fill the diabatic state populations
    for (int i = 0; i < opt.nstate; i++) {
        Pd.col(1 + i * (opt.nstate + 1)) = sdt(Eigen::all, Eigen::VectorXd::LinSpaced(opt.trajs, 0, opt.trajs - 1) * opt.nstate + Eigen::VectorXd::Constant(opt.trajs, i)).rowwise().sum() / opt.trajs;
    }

    // fill the adiabatic state populations
    for (int i = 0; i < opt.nstate; i++) {
        Pa.col(1 + i * (opt.nstate + 1)) = sat(Eigen::all, Eigen::VectorXd::LinSpaced(opt.trajs, 0, opt.trajs - 1) * opt.nstate + Eigen::VectorXd::Constant(opt.trajs, i)).rowwise().sum() / opt.trajs;
    }

    // print the new line
    std::cout << std::endl;

    // print the diabatic populations
    for (int i = 0; i < opt.nstate; i++) {
        std::printf("DIABATIC STATE %d POP: %.14f\n", i, Pd.bottomRows(1).rightCols(opt.nstate * opt.nstate).reshaped(opt.nstate, opt.nstate)(i, i));
    } std::cout << std::endl;

    // print the adiabatic populations
    for (int i = 0; i < opt.nstate; i++) {
        std::printf("ADIABATIC STATE %d POP: %.14f\n", i, Pa.bottomRows(1).rightCols(opt.nstate * opt.nstate).reshaped(opt.nstate, opt.nstate)(i, i));
    } std::cout << std::endl;

    // start the timer for writing the populations
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of the population writing
    std::cout << "POPULATION & HOPPING GEOMETRIES WRITING: " << std::flush;
        
    // write the populations to disk
    torch::WriteMatrix(std::string("P_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_DIA.mat", Pd); torch::WriteMatrix(std::string("P_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_ADIA.mat", Pa);

    // write the hopping geometries to disk
    torch::WriteMatrix(std::string("HG_LZ-") + (opt.adiabatic ? "ADIA" : "DIA") + "_.mat", Eigen::Map<Eigen::MatrixXd>(hg.data(), hg.size(), 1));

    // print the time for writing the populations
    std::cout << eltime(timers.at(1)) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << eltime(timers.at(0)) << std::endl;
}

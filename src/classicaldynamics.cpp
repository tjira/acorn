#include "classicaldynamics.h"

Eigen::MatrixXd ClassicalDynamics::evaluate_potential(std::vector<Expression>& potential_expressions, const Eigen::VectorXd& position) const {
    // deine the potential matrix
    Eigen::MatrixXd potential(input.potential.size(), input.potential.size());

    // evaluate the potential matrix
    for (size_t i = 0; i < input.potential.size(); i++) {
        for (size_t j = 0; j < input.potential.size(); j++) {
            potential(i, j) = potential_expressions.at(i * input.potential.size() + j).evaluate(position)(0);
        }
    }

    // return the potential matrix
    return potential;
}

Eigen::MatrixXd ClassicalDynamics::evaluate_potential_derivative(std::vector<Expression>& potential_expressions, const Eigen::VectorXd& position) const {
    // calculate the potential at some small spatial offset
    Eigen::MatrixXd potential_minus = evaluate_potential(potential_expressions, position.array() - 0.0001), potential_plus = evaluate_potential(potential_expressions, position.array() + 0.0001);

    // diagonalize the potentials if adiabatic dynamics requested
    if (input.adiabatic) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_minus(potential_minus), solver_plus(potential_plus); potential_minus = solver_minus.eigenvalues().asDiagonal(), potential_plus = solver_plus.eigenvalues().asDiagonal();
    }

    // return the potential derivative
    return (potential_plus - potential_minus) / 0.0002;
}

std::vector<Expression> ClassicalDynamics::get_potential_expression(const Wavefunction& initial_diabatic_wavefunction) const {
    // define the potential energy expressions
    std::vector<Expression> potential_expressions; potential_expressions.reserve(input.potential.size() * input.potential.size());

    // fill the potential energy expressions
    for (size_t i = 0; i < input.potential.size(); i++) {
        for (size_t j = 0; j < input.potential.size(); j++) {
            potential_expressions.emplace_back(input.potential.at(i).at(j), initial_diabatic_wavefunction.get_variables());
        }
    }

    // return the potential expressions
    return potential_expressions;
}

void ClassicalDynamics::print_iteration(int trajectory, int iteration, const std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXi>& trajectory_data, double mass) const {
    // extract the data
    const auto& [potential, position, velocity, states] = trajectory_data;

    // calculate the potential and kinetic energy
    double potential_energy = potential(states(iteration), states(iteration)), kinetic_energy = 0.5 * mass * velocity.squaredNorm();

    // print the fixed size information
    std::printf("%8d %8d %20.14f %20.14f %20.14f %5d [", trajectory, iteration, potential_energy, kinetic_energy, potential_energy + kinetic_energy, states(iteration) + 1);

    // print the position
    for (int k = 0; k < position.cols(); k++) {std::printf("%s%8.3f", k ? ", " : "", position(k));} std::printf("] [");

    // print the momentum
    for (int k = 0; k < velocity.cols(); k++) {std::printf("%s%8.3f", k ? ", " : "", velocity(k) * mass);} std::printf("]\n");
}

void ClassicalDynamics::run(const Input::Wavefunction& initial_diabatic_wavefunction_input) const {
    // define the initial wavefunction and extract the fourier grid
    auto [initial_diabatic_wavefunction, grid] = Wavefunction::Initialize(initial_diabatic_wavefunction_input); Eigen::MatrixXd fourier_grid = initial_diabatic_wavefunction.get_fourier_grid();

    // define the potential expressions and extract mass
    std::vector<Expression> potential_expressions = get_potential_expression(initial_diabatic_wavefunction); double mass = initial_diabatic_wavefunction.get_mass();

    // extract the initial position and momentum
    Eigen::VectorXd initial_position = initial_diabatic_wavefunction.position(grid), initial_momentum = initial_diabatic_wavefunction.momentum(fourier_grid);

    // define the initial state and the trajectory data
    int initial_diabatic_state, initial_adiabatic_state; std::vector<TrajectoryData> trajectory_data_vector(input.trajectories);

    // store the initial diabatic state
    initial_diabatic_wavefunction.get_density().diagonal().maxCoeff(&initial_diabatic_state);

    // define the diagonalizer for the initial potential
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> initial_potential_solver(evaluate_potential(potential_expressions, initial_position));

    // extract the initial adiabatic state
    (initial_potential_solver.eigenvectors() * initial_diabatic_wavefunction.get_density() * initial_potential_solver.eigenvectors().transpose()).diagonal().maxCoeff(&initial_adiabatic_state);

    // print the header with fixed length
    std::printf("\nCLASSICAL DYNAMICS\n%8s %8s %20s %20s %20s %5s", "TRAJ", "ITER", "EPOT", "EKIN", "ETOT", "STATE");

    // print the variable length header
    std::printf(" %*s %*s\n", (int)grid.cols() * 10, "POSITION", (int)grid.cols() * 10, "MOMENTUM");

    // loop over every trajectory
    for (int i = 0; i < input.trajectories; i++) {

        // distributions for position and momentum along with an uniform unitary distribution
        std::uniform_real_distribution<double> dist(0, 1); std::normal_distribution<double> positiondist(initial_position(0), 0.5), momentumdist(initial_momentum(0), 1.0);

        // random number generator and state vector with initial condition
        std::mt19937 mt(input.trajectories * (input.seed + i) + 1); Eigen::VectorXi state(input.iterations + 1); state(0) = input.adiabatic ? initial_adiabatic_state : initial_diabatic_state;

        // define and initialize the hopping algorithms
        LandauZener landauzener(input.adiabatic);

        // define the initial conditions
        Eigen::MatrixXd position(input.iterations + 1, grid.cols()), velocity(input.iterations + 1, grid.cols()), acceleration(input.iterations + 1, grid.cols());

        // fill the initial conditions
        position.row(0).fill(positiondist(mt)), velocity.row(0).fill(momentumdist(mt) / mass), acceleration.row(0).fill(0);

        // define a vector to contain diabatic potentials at each point
        std::vector<Eigen::MatrixXd> diabatic_potential_vector(input.iterations + 1), adiabatic_potential_vector(input.iterations + 1);

        // propagate the current trajectory
        for (int j = 0; j < input.iterations + 1; j++) {

            // evaluate the force
            Eigen::VectorXd force(grid.cols()); if (j) force(0) = -evaluate_potential_derivative(potential_expressions, position.row(j - 1))(state(j - 1), state(j - 1));

            // calculate the velocity and accceleration
            if (j) {acceleration.row(j) = force / mass; velocity.row(j) = velocity.row(j - 1) + 0.5 * (acceleration.row(j - 1) + acceleration.row(j)) * input.time_step;}

            // move the system and set the next state to the same as before
            if (j) {position.row(j) = position.row(j - 1) + input.time_step * (velocity.row(j) + 0.5 * acceleration.row(j) * input.time_step), state(j) = state(j - 1);}

            // calculate the potential at the current point and assign it to the container
            Eigen::MatrixXd potential = evaluate_potential(potential_expressions, position.row(j)); diabatic_potential_vector.at(j) = potential;

            // adiabatize the potential if requested
            if (input.adiabatic) {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(potential); potential = solver.eigenvalues().asDiagonal(), adiabatic_potential_vector.at(j) = potential;}

            // get the new state if we are at the crossing point according to the Landau-Zener algorithm
            int new_state = state(j); if (j > 1) new_state = landauzener.jump(input.adiabatic ? adiabatic_potential_vector : diabatic_potential_vector, j, state(j), input.time_step, dist(mt));

            // update the velocity and the state
            if (double new_velocity = velocity.row(j).squaredNorm() - 2 * (potential(new_state, new_state) - potential(state(j), state(j))) / mass; new_state != state(j) && new_velocity > 0) {
                velocity.row(j).array() = std::sqrt(velocity.row(j).squaredNorm() - 2 * (potential(new_state, new_state) - potential(state(j), state(j))) / mass); state(j) = new_state;
            }

            // print the iteration info
            if (j % input.log_interval == 0 || (j && state(j) != state(j - 1))) print_iteration(i + 1, j, {potential, position.row(j), velocity.row(j), state}, mass);
        }

        // set the trajectory data
        trajectory_data_vector.at(i) = {state, position, velocity, diabatic_potential_vector, adiabatic_potential_vector};
    }

    // create the vector of populations
    Eigen::VectorXd populations = Eigen::VectorXd::Zero(input.potential.size());

    // fill the vector with the populations
    for (int i = 0; i < input.trajectories; i++) {populations(trajectory_data_vector.at(i).state(input.iterations))++;} populations /= input.trajectories;

    // print the final populations
    for (size_t i = 0; i < input.potential.size(); i++) std::printf("%sFINAL %sDIABATIC POPULATION OF STATE %02lu: %.10f\n", i ? "" : "\n", input.adiabatic ? "A" : "", i + 1, populations(i));

    // start the timer nad print the export information
    Timepoint export_timer = Timer::Now(); std::printf("\nEXPORTING CLASSICAL DYNAMICS DATA: "); std::flush(std::cout);

    // export the trajectory data
    Export::ClassicalTrajectories(input, trajectory_data_vector, mass);

    // print the export time
    std::printf("%s\n", Timer::Format(Timer::Elapsed(export_timer)).c_str());
}

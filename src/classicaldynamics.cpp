#include "classicaldynamics.h"

Eigen::MatrixXd ClassicalDynamics::evaluate_potential(std::vector<std::vector<Expression>>& potential_expressions, const Eigen::VectorXd& position) const {
    // deine the potential matrix
    Eigen::MatrixXd potential(input.potential.size(), input.potential.size());

    // get the thread id
    #ifdef _OPENMP
    int id = omp_get_thread_num();
    #else
    int id = 0;
    #endif

    // evaluate the potential matrix
    for (size_t i = 0; i < input.potential.size(); i++) {
        for (size_t j = 0; j < input.potential.size(); j++) {
            potential(i, j) = potential_expressions.at(id).at(i * input.potential.size() + j).evaluate(position.transpose())(0);
        }
    }

    // return the potential matrix
    return potential;
}

Eigen::VectorXd ClassicalDynamics::evaluate_potential_derivative(std::vector<std::vector<Expression>>& potential_expressions, const Eigen::VectorXd& position, int state) const {
    // define the gradient vector
    Eigen::VectorXd gradient(position.size());

    // loop over every coordinate
    for (int i = 0; i < position.size(); i++) {

        // define the position with the small offset
        Eigen::VectorXd position_minus = position, position_plus = position; position_minus(i) -= 0.000001, position_plus(i) += 0.000001;

        // calculate the potential at the offset positions
        Eigen::MatrixXd potential_minus = evaluate_potential(potential_expressions, position_minus), potential_plus = evaluate_potential(potential_expressions, position_plus);

        // diagonalize the potentials if adiabatic dynamics requested
        if (input.adiabatic) {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_minus(potential_minus), solver_plus(potential_plus); potential_minus = solver_minus.eigenvalues().asDiagonal(), potential_plus = solver_plus.eigenvalues().asDiagonal();
        }

        // calculate the potential derivative
        gradient(i) = (potential_plus(state, state) - potential_minus(state, state)) / 0.000002;
    }

    // return the gradient
    return gradient;
}

std::vector<std::vector<Expression>> ClassicalDynamics::get_potential_expression(int dims) const {
    // define the potential energy expressions
    std::vector<std::vector<Expression>> potential_expressions(nthread); for (int i = 0; i < nthread; i++) potential_expressions.at(i).reserve(input.potential.size() * input.potential.size());

    // create the container for the variables
    std::vector<std::string> variables(dims);

    // fill the variables
    for (int i = 0; i < dims; i++) variables.at(i) = dims > 3 ? "r" + std::to_string(i + 1) : std::vector<std::string>{"x", "y", "z"}.at(i);

    // fill the potential energy expressions
    for (size_t i = 0; i < input.potential.size(); i++) {
        for (size_t j = 0; j < input.potential.size(); j++) {
            for (int k = 0; k < nthread; k++) potential_expressions.at(k).emplace_back(input.potential.at(i).at(j), variables, input.constants);
        }
    }

    // return the potential expressions
    return potential_expressions;
}

void ClassicalDynamics::print_iteration(int trajectory, int iteration, const std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXi>& trajectory_data, double mass) const {
    // extract the data and define the fixed size information
    const auto& [potential, position, velocity, states] = trajectory_data; char buffer[500];

    // calculate the potential and kinetic energy
    double potential_energy = potential(states(iteration), states(iteration)), kinetic_energy = 0.5 * mass * velocity.squaredNorm();

    // print the fixed size information
    std::sprintf(buffer, "%8d %8d %20.14f %20.14f %20.14f %5d [", trajectory, iteration, potential_energy, kinetic_energy, potential_energy + kinetic_energy, states(iteration) + 1);

    // print the position
    for (int k = 0; k < std::min((int)position.rows(), 3); k++) std::sprintf(buffer + strlen(buffer), "%s%8.3f", k ? ", " : "", position(k));

    // print the brackets
    std::sprintf(buffer + strlen(buffer), "%s] [" , position.rows() > 3 ? ", ..." : "");

    // print the momentum
    for (int k = 0; k < std::min((int)velocity.rows(), 3); k++) std::sprintf(buffer + strlen(buffer), "%s%8.3f", k ? ", " : "", velocity(k) * mass);

    // print the brackets
    std::sprintf(buffer + strlen(buffer), "%s]\n", velocity.rows() > 3 ? ", ..." : "");

    // print the buffer
    std::printf("%s", buffer);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd> ClassicalDynamics::read_initial_conditions(const Input::ClassicalDynamics::InitialConditions& ic) const {
    // throw an error if different lengths
    if (ic.positions.size() != ic.momenta.size() || ic.position_distribution.size() != ic.momentum_distribution.size()) throw std::runtime_error("THE NUMBER OF POSITIONS AND MOMENTA MUST BE THE SAME IN THE IC");

    // extract the dimension and initialize the random number generator
    int dims = ic.positions.empty() ? ic.position_distribution.size() : ic.positions.at(0).size(); std::mt19937 mt(input.seed); 

    // define the initial position, momentum and state vectors
    Eigen::MatrixXd initial_positions(input.trajectories, dims), initial_momenta(input.trajectories, dims); Eigen::VectorXd initial_diabatic_states(input.trajectories);

    // fill the initial conditions from the input
    for (int i = 0; i < input.trajectories && !ic.positions.empty(); i++) for (int j = 0; j < dims; j++) initial_positions(i, j)    = ic.positions.at(i).at(j);
    for (int i = 0; i < input.trajectories && !ic.positions.empty(); i++) for (int j = 0; j < dims; j++) initial_momenta(i, j)      =   ic.momenta.at(i).at(j);

    // fill the initial states
    std::fill(initial_diabatic_states.begin(), initial_diabatic_states.end(), ic.diabatic_state - 1);

    // return the initial conditions
    return {initial_positions, initial_momenta, initial_diabatic_states};
}

void ClassicalDynamics::run(const Input::Wavefunction& initial_diabatic_wavefunction_input) const {
    // check if the initial conditions are provided
    bool initial_conditions = (!input.initial_conditions.positions.empty() && !input.initial_conditions.momenta.empty()) || (!input.initial_conditions.position_distribution.empty() && !input.initial_conditions.momentum_distribution.empty());

    // throw an error if the specified initial state is adiabatic and diabatic dynamics is requested
    if (initial_conditions && !input.adiabatic && input.initial_conditions.adiabatic_state != -1) throw std::runtime_error("THE INITIAL ADIABATIC STATE CANNOT BE SPECIFIED IN DIABATIC DYNAMICS");

    // define the initial diabatic wavefunction, grids and initial confition variables
    Eigen::VectorXd initial_position, initial_momentum; int initial_diabatic_state = input.initial_conditions.diabatic_state - 1; double mass = input.initial_conditions.mass;

    // define the initial position, momentum and diabatic state for all trajectories
    Eigen::MatrixXd initial_positions, initial_momenta; Eigen::VectorXd initial_diabatic_states;

    // if no IC provided
    if (!initial_conditions) {

        // initialize the initial wavefunction and fourier grid
        auto [initial_diabatic_wavefunction, grid] = Wavefunction::Initialize(initial_diabatic_wavefunction_input); Eigen::MatrixXd fourier_grid = initial_diabatic_wavefunction.get_fourier_grid();

        // extract the initial position and momentum
        initial_position = initial_diabatic_wavefunction.position(grid), initial_momentum = initial_diabatic_wavefunction.momentum(fourier_grid);

        // extract the initial diabatic state and the fourier grid
        initial_diabatic_wavefunction.get_density().diagonal().maxCoeff(&initial_diabatic_state), mass = initial_diabatic_wavefunction.get_mass();
    }

    // if IC provided
    if (initial_conditions && !input.initial_conditions.positions.empty()) std::tie(initial_positions, initial_momenta, initial_diabatic_states) = read_initial_conditions(input.initial_conditions);

    // extract the number of dimensions
    int dims = initial_conditions ? (input.initial_conditions.positions.empty() ? input.initial_conditions.position_distribution.size() : initial_positions.cols()) : initial_position.rows();

    // define the potential expressions and trajectory data vector
    std::vector<std::vector<Expression>> potential_expressions = get_potential_expression(dims); std::vector<TrajectoryData> trajectory_data_vector(input.trajectories);

    // print the header with fixed length
    std::printf("CLASSICAL DYNAMICS\n%8s %8s %20s %20s %20s %5s", "TRAJ", "ITER", "EPOT", "EKIN", "ETOT", "STATE");

    // print the variable length header
    std::printf(" %*s %*s\n", std::min(dims, 3) * 10 + (dims > 3 ? 5 : 0), "POSITION", std::min(dims, 3) * 10 + (dims > 3 ? 5 : 0), "MOMENTUM");

    // loop over every trajectory
    #pragma omp parallel for num_threads(nthread) shared(trajectory_data_vector, potential_expressions)
    for (int i = 0; i < input.trajectories; i++) {

        // define the trajectory seed and initial adiabatic seed
        int seed = input.trajectories * (input.seed + i) + 1, initial_adiabatic_state = input.initial_conditions.adiabatic_state - 1;

        // define the mersenne twister for the trajectory and obtain initial conditions
        std::mt19937 mt(seed); auto [position, velocity, acceleration] = sample_initial_conditions(input.initial_conditions, initial_position, initial_momentum, mt, mass);

        // assign the fixed initial conditions if they were provided
        if (initial_positions.size() > 0) position.row(0) = initial_positions.row(i), velocity.row(0) = initial_momenta.row(i) / mass;

        // define the diagonalizer for the initial potential
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> initial_potential_solver(evaluate_potential(potential_expressions, position.row(0)));

        // if the adiabatic state is not specified, find the one with the highest population
        if (!initial_conditions || (initial_conditions && input.initial_conditions.adiabatic_state < 0)) {

            // form the initial diabatic density matrix
            Eigen::MatrixXd initial_diabatic_density = Eigen::MatrixXd::Zero(input.potential.size(), input.potential.size()); initial_diabatic_density(initial_diabatic_state, initial_diabatic_state) = 1;

            // extract the initial adiabatic state
            (initial_potential_solver.eigenvectors() * initial_diabatic_density * initial_potential_solver.eigenvectors().transpose()).diagonal().maxCoeff(&initial_adiabatic_state);
        }

        // define the state vector with initial condition
        Eigen::VectorXi state(input.iterations + 1); state(0) = input.adiabatic ? initial_adiabatic_state : initial_diabatic_state;

        // define and initialize the surface hopping algorithms
        FewestSwitches fewestswitches(input.surface_hopping, input.adiabatic, seed); LandauZener landauzener(input.surface_hopping, input.adiabatic, seed);

        // define the current fssh population vector and TDC
        Eigen::VectorXcd fssh_population = Eigen::VectorXcd::Zero(input.potential.size()); fssh_population(state(0)) = 1;

        // define the vector of hopping geometries and times
        std::vector<Eigen::VectorXd> hopping_geometry_vector; std::vector<double> hopping_time_vector;

        // define a vector to contain diabatic potentials and eigenvectors phi
        std::vector<Eigen::MatrixXd> diabatic_potential_vector(input.iterations + 1), adiabatic_potential_vector(input.iterations + 1), phi_vector(input.iterations + 1), tdc_vector(input.iterations);

        // propagate the current trajectory
        for (int j = 0; j < input.iterations + 1; j++) {

            // evaluate the force
            Eigen::VectorXd force(position.cols()); if (j) force = -evaluate_potential_derivative(potential_expressions, position.row(j - 1), state(j - 1));

            // calculate the velocity and accceleration
            if (j) {acceleration.row(j) = force / mass; velocity.row(j) = velocity.row(j - 1) + 0.5 * (acceleration.row(j - 1) + acceleration.row(j)) * input.time_step;}

            // move the system and set the next state to the same as before
            if (j) {position.row(j) = position.row(j - 1) + input.time_step * (velocity.row(j) + 0.5 * acceleration.row(j) * input.time_step), state(j) = state(j - 1);}

            // calculate the potential at the current point and assign it to the container and create the new state variable
            Eigen::MatrixXd potential = evaluate_potential(potential_expressions, position.row(j)), tdc; diabatic_potential_vector.at(j) = potential; int new_state = state(j); 

            // adiabatization block
            if (input.adiabatic) {

                // diagonalize the potential and assign the eigenvectors to the phi vector and the eigenvalues to the potential
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(potential); potential = solver.eigenvalues().asDiagonal(), adiabatic_potential_vector.at(j) = potential; phi_vector.at(j) = solver.eigenvectors();

                // flip the eigenvectors if the diagonalization algorithm found the opposite one
                if (j) for (int k = 0; k < phi_vector.at(j).cols(); k++) if ((phi_vector.at(j - 1).col(k).transpose() * phi_vector.at(j).col(k))(0) < 0) phi_vector.at(j).col(k) *= -1;
            }

            // get the new state if we are at the crossing point according to the Landau-Zener algorithm
            if (input.surface_hopping.type == "landau-zener" && j > 1) {
                new_state = landauzener.jump(input.adiabatic ? adiabatic_potential_vector : diabatic_potential_vector, j, state(j), input.time_step);
            } else if (input.surface_hopping.type == "fewest-switches" && j) {
                std::tie(tdc, fssh_population, new_state) = fewestswitches.jump(fssh_population, phi_vector, adiabatic_potential_vector, j, state(j), input.time_step);
            }

            // update the velocity and the state
            if (double Ekin = 0.5 * mass * velocity.row(j).squaredNorm(); new_state != state(j) && Ekin - (potential(new_state, new_state) - potential(state(j), state(j))) > 0) {
                velocity.row(j).array() *= std::sqrt((Ekin - (potential(new_state, new_state) - potential(state(j), state(j)))) / Ekin); state(j) = new_state;
            }

            // append the hopping geometry and time
            if (j && state(j) != state(j - 1) && input.data_export.hopping_geometry) hopping_geometry_vector.push_back(position.row(j));
            if (j && state(j) != state(j - 1) && input.data_export.hopping_time    ) hopping_time_vector.push_back(j * input.time_step);

            // append the TDC
            if (j && (input.data_export.tdc || input.data_export.tdc_mean)) tdc_vector.at(j - 1) = tdc;
            
            // print the iteration info
            if ((j % input.log_interval_step == 0 || (j && state(j) != state(j - 1))) && (i ? i + 1 : i) % input.log_interval_traj == 0) {
                print_iteration(i + 1, j, {potential, position.row(j), velocity.row(j), state}, mass);
            }
        }

        // set the trajectory data
        trajectory_data_vector.at(i) = {state, position, velocity, diabatic_potential_vector, adiabatic_potential_vector, tdc_vector, hopping_geometry_vector, hopping_time_vector};
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

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> ClassicalDynamics::sample_initial_conditions(const Input::ClassicalDynamics::InitialConditions& ic, const Eigen::VectorXd& initial_position, const Eigen::VectorXd& initial_momentum, std::mt19937& mt, double mass) const {
    // define the output matrices
    Eigen::MatrixXd position, velocity, acceleration;

    // if the sampling is needed
    if (initial_position.size() > 0 || !ic.position_distribution.empty()) {

        // define the vector of distributions
        std::vector<std::normal_distribution<double>> position_dist(!ic.position_distribution.empty() ? ic.position_distribution.size() : initial_position.size());
        std::vector<std::normal_distribution<double>> momentum_dist(!ic.momentum_distribution.empty() ? ic.momentum_distribution.size() : initial_momentum.size());

        // fill the distributions for the initial wavefuction
        for (size_t i = 0; i < position_dist.size() && ic.position_distribution.empty(); i++) position_dist.at(i) = std::normal_distribution<double>(initial_position(i), 0.5);
        for (size_t i = 0; i < momentum_dist.size() && ic.momentum_distribution.empty(); i++) momentum_dist.at(i) = std::normal_distribution<double>(initial_momentum(i), 1.0);

        // fill the distributions for the initial conditions
        for (size_t i = 0; i < position_dist.size() && !ic.position_distribution.empty(); i++) position_dist.at(i) = std::normal_distribution<double>(ic.position_distribution.at(i).at(0), ic.position_distribution.at(i).at(1));
        for (size_t i = 0; i < momentum_dist.size() && !ic.momentum_distribution.empty(); i++) momentum_dist.at(i) = std::normal_distribution<double>(ic.momentum_distribution.at(i).at(0), ic.momentum_distribution.at(i).at(1));

        // define the initial conditions
        position.resize(input.iterations + 1, position_dist.size()), velocity.resize(input.iterations + 1, momentum_dist.size()), acceleration.resize(input.iterations + 1, momentum_dist.size());

        // set all to zero
        position.setZero(), velocity.setZero(), acceleration.setZero();

        // fill the initial conditions for position and momentum
        for (size_t i = 0; i < position_dist.size(); i++) position(0, i) = position_dist.at(i)(mt);
        for (size_t i = 0; i < momentum_dist.size(); i++) velocity(0, i) = momentum_dist.at(i)(mt);

        // momentum to velocity
        velocity.row(0) /= mass;

    // if the sampling is not nesessary and explicit initial conditions are provided, set the correct shape of the matrices
    } else position = Eigen::MatrixXd::Zero(input.iterations + 1, ic.positions.at(0).size()), velocity = position, acceleration = position;

    // return the initial conditions
    return {position, velocity, acceleration};
}

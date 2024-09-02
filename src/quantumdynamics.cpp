#include "quantumdynamics.h"

Eigen::MatrixXd QuantumDynamics::get_diabatic_potential(const Eigen::MatrixXd& grid, const std::vector<std::string>& variables) const {
    // create the container for the diabatic potential
    Eigen::MatrixXd diabatic_potential(grid.rows(), (int)std::pow(input.potential.size(), 2));

    // fill the diabatic potential
    for (size_t i = 0; i < input.potential.size(); i++) for (size_t j = 0; j < input.potential.size(); j++) {
        diabatic_potential.col(i * input.potential.size() + j) = Expression(input.potential.at(i).at(j), variables).evaluate(grid);
    }

    // return the diabatic potential
    return diabatic_potential;
}

std::vector<Eigen::MatrixXd> QuantumDynamics::get_transformation_matrices(const Eigen::MatrixXd& diabatic_potential) const {
    // define the vector of transformation matrices
    std::vector<Eigen::MatrixXd> transformation_matrices(diabatic_potential.rows(), Eigen::MatrixXd::Identity(input.potential.size(), input.potential.size()));

    // diagonalize the potential at each point
    for (int i = 0; i < diabatic_potential.rows(); i++) {

        // initialize the diabatic potential at the current coordinates
        Eigen::MatrixXd coordinate_potential = diabatic_potential.row(i).reshaped(input.potential.size(), input.potential.size());

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(coordinate_potential); Eigen::MatrixXd C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(input.potential.size(), 1);

        // calculate the overlap of eigenvectors
        for (size_t j = 0; j < input.potential.size(); j++) {
            overlap(j) = transformation_matrices.at(std::max(i - 1, 0)).col(j).transpose() * C.col(j); overlap(j) = overlap(j) < 0 ? -1 : 1;
        }

        // maximize the overlap with the previuos transformation
        transformation_matrices.at(i) = C * overlap.asDiagonal();
    }

    // return the transformation matrices
    return transformation_matrices;
}

void QuantumDynamics::run(const Wavefunction& initial_diabatic_wavefunction) const {
    // generate the independent variable grid and the fourier grid
    Eigen::MatrixXd grid = initial_diabatic_wavefunction.get_grid(), fourier_grid = initial_diabatic_wavefunction.get_fourier_grid();

    // evaluate the potential on the precalculated grid
    Eigen::MatrixXd diabatic_potential = get_diabatic_potential(grid, initial_diabatic_wavefunction.get_variables());

    // calculate the adiabatic transformation matrices
    std::vector<Eigen::MatrixXd> transformation_matrices = get_transformation_matrices(diabatic_potential);

    // calculate the propagators for the imaginary propagation
    auto [real_propagators, fourier_propagators] = initial_diabatic_wavefunction.propagators(diabatic_potential, fourier_grid, std::complex<double>(1, 0), input.time_step);

    // create a vector of each state that will get propagated imaginary times
    std::vector<Wavefunction> imaginary_diabatic_states(input.imaginary, initial_diabatic_wavefunction);

    // loop over all imaginary propagated states
    for (int i = 0; i < input.imaginary; i++) {

        // calculate the initial energy and define the container for the diabatic wavefunction trajectory
        double energy = imaginary_diabatic_states.at(i).energy(diabatic_potential, fourier_grid); std::vector<Wavefunction> diabatic_wavefunction_trajectory(input.iterations + 1);

        // assign the zeroth interation wavefunction to the trajectory container
        if (input.export_data) diabatic_wavefunction_trajectory.at(0) = imaginary_diabatic_states.at(i);

        // print the header and the zeroth iteration
        std::printf("\nIMAGINARY PROPAGATION OF STATE %d\n%6s %20s %8s %12s\n%6d %20.14f %.2e 00:00:00.000\n", i + 1, "ITER", "ENERGY", "|dE|", "TIME", 0, energy, 0.0);

        // perform the propagation in imaginary time
        for (int j = 0; j < input.iterations; j++) {

            // start the timer
            Timepoint iteration_timer = Timer::Now();

            // save the previous energy and propagate the state
            double energy_old = energy; imaginary_diabatic_states.at(i) = imaginary_diabatic_states.at(i).propagated(real_propagators, fourier_propagators);

            // orthogonalize the imaginary propagated state
            for (int k = 0; k < i; k++) {
                imaginary_diabatic_states.at(i) = imaginary_diabatic_states.at(i) - imaginary_diabatic_states.at(k) * imaginary_diabatic_states.at(k).overlap(imaginary_diabatic_states.at(i));
            }

            // normalize the wavefunction and calculate the new energy
            imaginary_diabatic_states.at(i) = imaginary_diabatic_states.at(i).get_normalized(), energy = imaginary_diabatic_states.at(i).energy(diabatic_potential, fourier_grid);

            // append the wavefunction to the container
            if (input.export_data) diabatic_wavefunction_trajectory.at(j + 1) = imaginary_diabatic_states.at(i);

            // print the iteration information
            std::printf("%6d %20.14f %.2e %s\n", j + 1, energy, std::abs(energy_old - energy), Timer::Format(Timer::Elapsed(iteration_timer)).c_str());
        }

        // export the wavefunction
        if (input.export_data) Export::WavefunctionTrajectory("PSI_DIA_IMAG_" + std::to_string(i + 1) + ".mat", diabatic_wavefunction_trajectory, grid);
    }

    // print the energies of the individual states
    for (int i = 0; i < input.imaginary; i++) std::printf("%sOPTIMIZED WAVEFUNCTION %02d ENERGY: %.14f\n", i ? "" : "\n", i + 1, imaginary_diabatic_states.at(i).energy(diabatic_potential, fourier_grid));

    // calculate the propagators for the real propagation
    std::tie(real_propagators, fourier_propagators) = initial_diabatic_wavefunction.propagators(diabatic_potential, fourier_grid, std::complex<double>(0, 1), input.time_step);

    // create a vector of each state that will get propagated real times
    std::vector<Wavefunction> real_diabatic_states(input.real, initial_diabatic_wavefunction);

    // assign the initial conditions from the imaginary propagation
    for (int i = 0; i < std::min(input.imaginary, input.real); i++) real_diabatic_states.at(i) = imaginary_diabatic_states.at(i);

    // loop over all real propagated states
    for (int i = 0; i < input.real; i++) {

        // calculate the initial energy and define the container for the diabatic wavefunction trajectory
        double energy = real_diabatic_states.at(i).energy(diabatic_potential, fourier_grid); std::vector<Wavefunction> diabatic_wavefunction_trajectory(input.iterations + 1);

        // assign the zeroth interation wavefunction to the trajectory container
        if (input.export_data) diabatic_wavefunction_trajectory.at(0) = real_diabatic_states.at(i);

        // define the density matrices in diabatic and adiabatic bases
        Eigen::MatrixXd density_diabatic(input.iterations + 1, (int)std::pow(input.potential.size(), 2)), density_adiabatic(input.iterations + 1, (int)std::pow(input.potential.size(), 2));

        // fill the initial diabatic and adiabatic density matrix
        density_diabatic.row(0) = real_diabatic_states.at(i).get_density().reshaped(), density_adiabatic.row(0) = real_diabatic_states.at(i).adiabatized(transformation_matrices).get_density().reshaped();

        // print the header and the zeroth iteration
        std::printf("\nREAL PROPAGATION OF STATE %d\n%6s %20s %8s %12s\n%6d %20.14f %.2e 00:00:00.000\n", i + 1, "ITER", "ENERGY", "|dE|", "TIME", 0, energy, 0.0);

        // perform the propagation in imaginary time
        for (int j = 0; j < input.iterations; j++) {

            // start the timer
            Timepoint iteration_timer = Timer::Now();

            // save the previous energy and propagate the state
            double energy_old = energy; real_diabatic_states.at(i) = real_diabatic_states.at(i).propagated(real_propagators, fourier_propagators);

            // calculate the new energy and append the wavefunction to the container
            energy = real_diabatic_states.at(i).energy(diabatic_potential, fourier_grid); if (input.export_data) diabatic_wavefunction_trajectory.at(j + 1) = real_diabatic_states.at(i);

            // calculate and assign the diabatic and adiabatic density matrix
            density_diabatic .row(j + 1) = real_diabatic_states.at(i)                                     .get_density().reshaped();
            density_adiabatic.row(j + 1) = real_diabatic_states.at(i).adiabatized(transformation_matrices).get_density().reshaped();

            // print the iteration information
            std::printf("%6d %20.14f %.2e %s\n", j + 1, energy, std::abs(energy_old - energy), Timer::Format(Timer::Elapsed(iteration_timer)).c_str());
        }

        // print the wavefunction diabatic state populations
        for (size_t j = 0; j < input.potential.size(); j++) {
            std::printf("%sFINAL DIABATIC STATE %02lu POPULATION: %.14f\n", j ? "" : "\n", j + 1, density_diabatic.bottomRows(1).reshaped(input.potential.size(), input.potential.size())(j, j));
        }

        // print the wavefunction adiabatic state populations
        for (size_t j = 0; j < input.potential.size(); j++) {
            std::printf("%sFINAL ADIABATIC STATE %02lu POPULATION: %.14f\n", j ? "" : "\n", j + 1, density_adiabatic.bottomRows(1).reshaped(input.potential.size(), input.potential.size())(j, j));
        }

        // export the density matrices and wavefunctions
        if (auto time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.iterations * input.time_step); input.export_data) {
            Export::EigenMatrixDouble     ("P_DIA_"        + std::to_string(i + 1) + ".mat", density_diabatic,                 time_variable);
            Export::EigenMatrixDouble     ("P_ADIA_"       + std::to_string(i + 1) + ".mat", density_adiabatic,                time_variable);
            Export::WavefunctionTrajectory("PSI_DIA_REAL_" + std::to_string(i + 1) + ".mat", diabatic_wavefunction_trajectory, grid         );
        }
    }

    // export the potential
    if (input.export_data) Export::EigenMatrixDouble("U_DIA.mat", diabatic_potential, grid);

    // print the end energies of the individual states
    for (int i = 0; i < input.real; i++) std::printf("%sFINAL WAVEFUNCTION %02d ENERGY: %.14f\n", i ? "" : "\n", i + 1, real_diabatic_states.at(i).energy(diabatic_potential, fourier_grid));
}

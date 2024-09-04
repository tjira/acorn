#include "quantumdynamics.h"

void QuantumDynamics::export_trajectory(int iteration, const std::vector<IterationData>& iteration_data, const Eigen::MatrixXd& grid, bool imaginary) const {
    // create the data matrix
    Eigen::MatrixXd data_matrix;

    // export the diabatic wavefunctions
    if (input.data_export.diabatic_wavefunction) {

        // create the container for the wavefunction trajectory
        data_matrix = Eigen::MatrixXd(iteration_data.front().diabatic_wavefunction.get_data().rows(), 2 * iteration_data.size() * iteration_data.front().diabatic_wavefunction.get_data().cols());

        // fill the data matrix with diabatic wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) {

            // extract the wavefunction data
            const Eigen::MatrixXcd& wavefunction_data = iteration_data.at(i).diabatic_wavefunction.get_data();

            // for each column of the data
            for (int j = 0; j < wavefunction_data.cols(); j++) {

                // extract the real and imaginary parts of the wavefunction column
                const Eigen::VectorXd& real_part = wavefunction_data.col(j).real();
                const Eigen::VectorXd& imag_part = wavefunction_data.col(j).imag();

                // copy the real and imaginary parts to the wavefunction matrix
                data_matrix.block(0, 2 * (i * wavefunction_data.cols() + j) + 0, wavefunction_data.rows(), 1) = real_part;
                data_matrix.block(0, 2 * (i * wavefunction_data.cols() + j) + 1, wavefunction_data.rows(), 1) = imag_part;
            }
        }

        // save the diabatic wavefunction
        Export::EigenMatrixDouble(std::string("PSI_DIABATIC_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, grid);
    }

    // export the adiabatic wavefunctions
    if (input.data_export.adiabatic_wavefunction) {

        // create the container for the wavefunction trajectory
        data_matrix = Eigen::MatrixXd(iteration_data.front().adiabatic_wavefunction.get_data().rows(), 2 * iteration_data.size() * iteration_data.front().adiabatic_wavefunction.get_data().cols());

        // fill the data matrix with diabatic wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) {

            // extract the wavefunction data
            const Eigen::MatrixXcd& wavefunction_data = iteration_data.at(i).adiabatic_wavefunction.get_data();

            // for each column of the data
            for (int j = 0; j < wavefunction_data.cols(); j++) {

                // extract the real and imaginary parts of the wavefunction column
                const Eigen::VectorXd& real_part = wavefunction_data.col(j).real();
                const Eigen::VectorXd& imag_part = wavefunction_data.col(j).imag();

                // copy the real and imaginary parts to the wavefunction matrix
                data_matrix.block(0, 2 * (i * wavefunction_data.cols() + j) + 0, wavefunction_data.rows(), 1) = real_part;
                data_matrix.block(0, 2 * (i * wavefunction_data.cols() + j) + 1, wavefunction_data.rows(), 1) = imag_part;
            }
        }

        // save the diabatic wavefunction
        Export::EigenMatrixDouble(std::string("PSI_ADIABATIC_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, grid);
    }

    // export the diabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(iteration_data.size(), 0, input.time_step * iteration_data.size()); input.data_export.diabatic_density) {

        // create the container for the density matrix
        data_matrix = Eigen::MatrixXd(iteration_data.size(), input.potential.size() * input.potential.size());

        // fill the data matrix with diabatic density
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).density_diabatic.reshaped();

        // save the diabatic wavefunction
        Export::EigenMatrixDouble(std::string("DENSITY_DIABATIC_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, time_variable);
    }

    // export the adiabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(iteration_data.size(), 0, input.time_step * iteration_data.size()); input.data_export.adiabatic_density) {

        // create the container for the density matrix
        data_matrix = Eigen::MatrixXd(iteration_data.size(), input.potential.size() * input.potential.size());

        // fill the data matrix with diabatic density
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).density_adiabatic.reshaped();

        // save the diabatic wavefunction
        Export::EigenMatrixDouble(std::string("DENSITY_ADIABATIC_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, time_variable);
    }

    // export the position
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(iteration_data.size(), 0, input.time_step * iteration_data.size()); input.data_export.position) {

        // create the container for the position
        data_matrix = Eigen::MatrixXd(iteration_data.size(), iteration_data.front().position.size());

        // fill the data matrix with position of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).position;

        // save the wavefunction position
        Export::EigenMatrixDouble(std::string("POSITION_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, time_variable);
    }

    // export the momentum
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(iteration_data.size(), 0, input.time_step * iteration_data.size()); input.data_export.momentum) {

        // create the container for the momentum
        data_matrix = Eigen::MatrixXd(iteration_data.size(), iteration_data.front().momentum.size());

        // fill the data matrix with momentum of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).momentum;

        // save the wavefunction momentum
        Export::EigenMatrixDouble(std::string("MOMENTUM_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, time_variable);
    }

    // export the energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(iteration_data.size(), 0, input.time_step * iteration_data.size()); input.data_export.energy) {

        // create the container for the energy
        data_matrix = Eigen::MatrixXd(iteration_data.size(), 1);

        // fill the data matrix with energy of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix(i) = iteration_data.at(i).energy;

        // save the wavefunction energy
        Export::EigenMatrixDouble(std::string("ENERGY_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, time_variable);
    }

    // export the autocorrelation function
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(iteration_data.size(), 0, input.time_step * iteration_data.size()); input.data_export.acf) {

        // create the container for the autocorrelation function
        data_matrix = Eigen::MatrixXd(iteration_data.size(), 2);

        // fill the data matrix with autocorrelation function of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix(i, 0) = iteration_data.at(i).acf.real(), data_matrix(i, 1) = iteration_data.at(i).acf.imag();

        // save the wavefunction autocorrelation function
        Export::EigenMatrixDouble(std::string("ACF_") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(iteration + 1) + ".mat", data_matrix, time_variable);
    }
}

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

std::tuple<std::vector<Eigen::MatrixXd>, Eigen::MatrixXd> QuantumDynamics::get_transformation_matrices(const Eigen::MatrixXd& diabatic_potential) const {
    // define the vector of transformation matrices
    std::vector<Eigen::MatrixXd> transformation_matrices(diabatic_potential.rows(), Eigen::MatrixXd::Identity(input.potential.size(), input.potential.size()));

    // define the adiabatic potential matrix
    Eigen::MatrixXd adiabatic_potential(diabatic_potential.rows(), input.potential.size());

    // diagonalize the potential at each point
    for (int i = 0; i < diabatic_potential.rows(); i++) {

        // initialize the diabatic potential at the current coordinates
        Eigen::MatrixXd coordinate_potential = diabatic_potential.row(i).reshaped(input.potential.size(), input.potential.size());

        // solve the eigenvalue problem and define overlap
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(coordinate_potential); Eigen::MatrixXd C = solver.eigenvectors(), eps = solver.eigenvalues(), overlap(input.potential.size(), 1);

        // addign the eigenvalues to the adiabatic potential
        adiabatic_potential.row(i) = eps.transpose();

        // calculate the overlap of eigenvectors
        for (size_t j = 0; j < input.potential.size(); j++) {
            overlap(j) = transformation_matrices.at(std::max(i - 1, 0)).col(j).transpose() * C.col(j); overlap(j) = overlap(j) < 0 ? -1 : 1;
        }

        // maximize the overlap with the previuos transformation
        transformation_matrices.at(i) = C * overlap.asDiagonal();
    }

    // return the transformation matrices
    return {transformation_matrices, adiabatic_potential};
}

void QuantumDynamics::print_iteration(int iteration, const IterationData& iteration_data, long elapsed) const {
    // print the main iteration info
    std::printf("%6d %20.14f %20.14f %20.14f %.2e %s [", iteration, iteration_data.energy, iteration_data.acf.real(), iteration_data.acf.imag(), iteration_data.energy_error, Timer::Format(elapsed).c_str());

    // print the wavefunction position
    for (int i = 0; i < iteration_data.position.size(); i++) {std::printf("%s%8.3f", i ? ", " : "", iteration_data.position(i));} std::printf("] [");

    // print the wavefunction momentum
    for (int i = 0; i < iteration_data.momentum.size(); i++) {std::printf("%s%8.3f", i ? ", " : "", iteration_data.momentum(i));} std::printf("] [");

    // print the diabaic density matrix
    for (size_t i = 0; i < input.potential.size(); i++) {std::printf("%s%8.3f", i ? ", " : "", iteration_data.density_diabatic(i, i));} std::printf("]\n");
}

void QuantumDynamics::run(const Input::Wavefunction& initial_diabatic_wavefunction_input) const {
    // initialize the timer for grid generation and print the grid generation timer
    Timepoint grid_generation_timer = Timer::Now(); std::printf("\nINITIAL GRIDS AND POTENTIALS GENERATION: "); std::flush(std::cout);

    // initialize the initial wavefunction and the grids
    auto [initial_diabatic_wavefunction, grid] = Wavefunction::Initialize(initial_diabatic_wavefunction_input); Eigen::MatrixXd fourier_grid = initial_diabatic_wavefunction.get_fourier_grid();

    // evaluate the potential on the precalculated grid
    Eigen::MatrixXd diabatic_potential = get_diabatic_potential(grid, initial_diabatic_wavefunction.get_variables());

    // calculate the adiabatic transformation matrices and potential
    auto [transformation_matrices, adiabatic_potential] = get_transformation_matrices(diabatic_potential);

    // export the potentials
    if (input.data_export.adiabatic_potential) Export::EigenMatrixDouble("POTENTIAL_ADIABATIC.mat", adiabatic_potential, grid);
    if (input.data_export.diabatic_potential ) Export::EigenMatrixDouble("POTENTIAL_DIABATIC.mat",  diabatic_potential,  grid);

    // print the grid generation timing
    std::printf("%s\n", Timer::Format(Timer::Elapsed(grid_generation_timer)).c_str());

    // create a vector of each state that will get propagated in imaginary time
    std::vector<Wavefunction> imaginary_diabatic_states(input.imaginary, initial_diabatic_wavefunction);

    // imaginary and real propagation
    for (int i = 0; i < 2; i++) {

        // define the state count, the imaginary flag and vector of energies
        int state_count = i ? input.real : input.imaginary; bool imaginary = !i; std::vector<double> energies(state_count); if (!state_count) continue;

        // initialize the timer for propagator generation and print the header
        Timepoint propagator_generation_timer = Timer::Now(); std::printf("%sREAL AND FOURIER PROPAGATORS GENERATION: ", !imaginary && input.imaginary ? "\n" : ""); std::flush(std::cout);

        // calculate the propagators
        auto [real_propagators, fourier_propagators] = initial_diabatic_wavefunction.propagators(diabatic_potential, fourier_grid, std::complex<double>(imaginary, !imaginary), input.time_step);

        // print the propagator generation timing
        std::printf("%s\n", Timer::Format(Timer::Elapsed(propagator_generation_timer)).c_str());

        // loop over all propagated states
        for (int j = 0; j < state_count; j++) {

            // define the iteration data structure
            std::vector<IterationData> iteration_data_vector(input.iterations + 1);

            // assign the initial conditions to the container of the wavefunction to propagate
            Wavefunction diabatic_wavefunction = input.imaginary > j ? imaginary_diabatic_states.at(j) : initial_diabatic_wavefunction, adiabatic_wavefunction;

            // print the header without the variable dimension variables
            std::printf("\n%s PROPAGATION OF STATE %d\n%6s %20s %20s %20s %8s %12s", imaginary ? "IMAGINARY" : "REAL", j + 1, "ITER", "ENERGY", "RE(ACF)", "IM(ACF)", "|dE|", "TIME");

            // print the variable length header
            std::printf(" %*s %*s %*s\n", (int)grid.cols() * 10, "POSITION", (int)grid.cols() * 10, "MOMENTUM", (int)(grid.cols() * input.potential.size()), "POPULATION");

            // perform the propagation in imaginary time
            for (int k = 0; k < input.iterations + 1; k++) {

                // start the timer, save the previous energy and define the iteration data
                Timepoint iteration_timer = Timer::Now(); double energy_old = k ? iteration_data_vector.at(k - 1).energy : 0; IterationData iteration_data;

                // propagate the wavefunction in diabatic basis
                if (k) diabatic_wavefunction = (input.data_export.diabatic_wavefunction ? iteration_data_vector.at(k - 1).diabatic_wavefunction : diabatic_wavefunction).propagated(real_propagators, fourier_propagators);

                // orthogonalize the imaginary propagated state
                for (int l = 0; l < j && imaginary; l++) diabatic_wavefunction = diabatic_wavefunction - imaginary_diabatic_states.at(l) * imaginary_diabatic_states.at(l).overlap(diabatic_wavefunction);

                // normalize the wavefunction if imaginary time propagation is performed
                if (imaginary) {diabatic_wavefunction = diabatic_wavefunction.get_normalized(); if (!k) imaginary_diabatic_states.at(j) = diabatic_wavefunction;}

                // transform the wavefunction to adiabatic basis
                if (input.data_export.adiabatic_wavefunction || input.data_export.adiabatic_density) adiabatic_wavefunction = diabatic_wavefunction.adiabatized(transformation_matrices);

                // calculate all the properties of the wavefunction
                iteration_data.energy   = diabatic_wavefunction.energy(diabatic_potential, fourier_grid), iteration_data.energy_error = std::abs(iteration_data.energy - energy_old);
                iteration_data.position = diabatic_wavefunction.position(grid),                           iteration_data.momentum     = diabatic_wavefunction.momentum(fourier_grid);

                // calculate the autocorrelation function
                iteration_data.acf = (input.imaginary > j ? imaginary_diabatic_states.at(j) : initial_diabatic_wavefunction).overlap(diabatic_wavefunction);

                // calculate the density matrices
                iteration_data.density_diabatic = diabatic_wavefunction.get_density(); if (input.data_export.adiabatic_density) iteration_data.density_adiabatic = adiabatic_wavefunction.get_density();

                // move the wavefunctions to the container if requested
                if (input.data_export.diabatic_wavefunction ) iteration_data.diabatic_wavefunction  = std::move(diabatic_wavefunction );
                if (input.data_export.adiabatic_wavefunction) iteration_data.adiabatic_wavefunction = std::move(adiabatic_wavefunction);

                // append the iteration data to the container and print the iteration
                iteration_data_vector.at(k) = iteration_data; print_iteration(k, iteration_data, Timer::Elapsed(iteration_timer));
            }

            // if the imaginary time propagation is performed, append the wafefunction to the optimized vector
            if (imaginary) imaginary_diabatic_states.at(j) = input.data_export.diabatic_wavefunction ? iteration_data_vector.back().diabatic_wavefunction : diabatic_wavefunction;

            // print the wavefunction diabatic state populations
            for (size_t k = 0; k < input.potential.size(); k++) {
                std::printf("%sFINAL DIABATIC  STATE %02lu POPULATION (%s): %.10f\n", k ? "" : "\n", k + 1, imaginary ? "ITP" : "RTP", iteration_data_vector.back().density_diabatic(k, k));
            }

            // start the timer nad print the export information
            Timepoint export_timer = Timer::Now(); std::printf("EXPORTING DATA FOR WAVEFUNCTION %03d (%s): ", j + 1, imaginary ? "ITP" : "RTP"); std::flush(std::cout);

            // export the wavefunction trajectory with data and save the energy
            export_trajectory(j, iteration_data_vector, grid, imaginary); energies.at(j) = iteration_data_vector.back().energy;

            // print the export time
            std::printf("%s\n", Timer::Format(Timer::Elapsed(export_timer)).c_str());
        }
        
        // print the final energies
        for (size_t j = 0; j < energies.size(); j++) {
            std::printf("%sFINAL ENERGY OF WAVEFUNCTION %02lu (%s): %20.14f\n", j ? "" : "\n", j + 1, imaginary ? "ITP" : "RTP", energies.at(j));
        }
    }
}

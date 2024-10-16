#include "classicaldynamics.h"
#include            "export.h"
#include   "quantumdynamics.h"

void Export::ClassicalTrajectories(const Input::ClassicalDynamics& input, const std::vector<ClassicalDynamics::TrajectoryData>& trajectory_data_vector, int mass) {
    // define the output data matrix and the eigenvalue solver
    Eigen::MatrixXd data_matrix; Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

    // get the surface hopping algorithm name
    std::string algorithm = input.surface_hopping.type == "landau-zener" ? "LZ" : input.surface_hopping.type == "fewest-switches" ? "FS" : "";

    // export the diabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.diabatic_population) {

        // create the matrix containing the diabatic populations
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, input.potential.size());

        // fill the matrix with the diabatic populations
        for (int i = 0; i < input.trajectories; i++) {
            for (int j = 0; j < input.iterations + 1; j++) {

                // define the pseudo density matrix
                Eigen::VectorXd pseudo_density = Eigen::VectorXd::Zero(input.potential.size()); pseudo_density(trajectory_data_vector.at(i).state(j)) = 1;

                // create the solver for the potential matrix
                if (input.adiabatic) solver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(trajectory_data_vector.at(i).diabatic_potential.at(j));

                // transform the pseudo density matrix
                if (input.adiabatic) pseudo_density = (solver.eigenvectors() * pseudo_density.asDiagonal() * solver.eigenvectors().transpose()).diagonal().transpose();

                // assign the adiabatic populations
                data_matrix.row(j) += pseudo_density;
            }
        }

        // export the diabatic populations
        Export::EigenMatrixDouble(std::string("POPULATION-DIABATIC_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the diabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.adiabatic_population) {

        // create the matrix containing the diabatic populations
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, input.potential.size());

        // fill the matrix with the diabatic populations
        for (int i = 0; i < input.trajectories; i++) {
            for (int j = 0; j < input.iterations + 1; j++) {

                // define the pseudo density matrix
                Eigen::VectorXd pseudo_density = Eigen::VectorXd::Zero(input.potential.size()); pseudo_density(trajectory_data_vector.at(i).state(j)) = 1;

                // create the solver for the potential matrix
                if (!input.adiabatic) solver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(trajectory_data_vector.at(i).diabatic_potential.at(j));

                // transform the pseudo density matrix
                if (!input.adiabatic) pseudo_density = (solver.eigenvectors().transpose() * pseudo_density.asDiagonal() * solver.eigenvectors()).diagonal().transpose();

                // assign the adiabatic populations
                data_matrix.row(j) += pseudo_density;
            }
        }

        // export the diabatic populations
        Export::EigenMatrixDouble(std::string("POPULATION-ADIABATIC_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the position
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.position) {

        // create the matrix containing the position
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, trajectory_data_vector.size());

        // fill the matrix with the position
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) data_matrix(j, i) = trajectory_data_vector.at(i).position.row(j)(0);

        // export the position
        Export::EigenMatrixDouble(std::string("POSITION_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix, time_variable);
    }

    // export the position mean
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.position_mean) {

        // create the matrix containing the position mean
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, 1);

        // fill the matrix with the position
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) data_matrix(j, 0) += trajectory_data_vector.at(i).position.row(j)(0);

        // export the position
        Export::EigenMatrixDouble(std::string("POSITION-MEAN_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the momentum
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.momentum) {

        // create the matrix containing the momentum
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, trajectory_data_vector.size());

        // fill the matrix with the momentum
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) data_matrix(j, i) = trajectory_data_vector.at(i).velocity.row(j)(0) * mass;

        // export the momentum
        Export::EigenMatrixDouble(std::string("MOMENTUM_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix, time_variable);
    }

    // export the momentum mean
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.momentum_mean) {

        // create the matrix containing the momentum
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, 1);

        // fill the matrix with the momentum
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) data_matrix(j, 0) += trajectory_data_vector.at(i).velocity.row(j)(0) * mass;

        // export the momentum
        Export::EigenMatrixDouble(std::string("MOMENTUM-MEAN_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the total energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.total_energy) {

        // create the matrix containing the energy
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, trajectory_data_vector.size());

        // fill the matrix with the energy
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) {

            // extract the potential matrix
            const Eigen::MatrixXd& potential = input.adiabatic ? trajectory_data_vector.at(i).adiabatic_potential.at(j) : trajectory_data_vector.at(i).diabatic_potential.at(j);

            // extract the state and assign the energy
            int state = trajectory_data_vector.at(i).state(j); data_matrix.row(j)(i) = potential(state, state) + 0.5 * mass * trajectory_data_vector.at(i).velocity.row(j).squaredNorm();
        }

        // export the energy
        Export::EigenMatrixDouble(std::string("TOTAL-ENERGY_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix, time_variable);
    }

    // export the total energy mean
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.total_energy_mean) {

        // create the matrix containing the energy
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, 1);

        // fill the matrix with the energy
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) {

            // extract the potential matrix
            const Eigen::MatrixXd& potential = input.adiabatic ? trajectory_data_vector.at(i).adiabatic_potential.at(j) : trajectory_data_vector.at(i).diabatic_potential.at(j);

            // extract the state and assign the energy
            int state = trajectory_data_vector.at(i).state(j); data_matrix(j, 0) += potential(state, state) + 0.5 * mass * trajectory_data_vector.at(i).velocity.row(j).squaredNorm();
        }

        // export the energy mean
        Export::EigenMatrixDouble(std::string("TOTAL-ENERGY-MEAN_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the potential energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.potential_energy) {

        // create the matrix containing the energy
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, trajectory_data_vector.size());

        // fill the matrix with the energy
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) {

            // extract the potential matrix
            const Eigen::MatrixXd& potential = input.adiabatic ? trajectory_data_vector.at(i).adiabatic_potential.at(j) : trajectory_data_vector.at(i).diabatic_potential.at(j);

            // extract the state and assign the energy
            int state = trajectory_data_vector.at(i).state(j); data_matrix.row(j)(i) = potential(state, state);
        }

        // export the energy
        Export::EigenMatrixDouble(std::string("POTENTIAL-ENERGY_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix, time_variable);
    }

    // export the potential energy mean
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.potential_energy_mean) {

        // create the matrix containing the energy
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, 1);

        // fill the matrix with the energy
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) {

            // extract the potential matrix
            const Eigen::MatrixXd& potential = input.adiabatic ? trajectory_data_vector.at(i).adiabatic_potential.at(j) : trajectory_data_vector.at(i).diabatic_potential.at(j);

            // extract the state and assign the energy
            int state = trajectory_data_vector.at(i).state(j); data_matrix(j, 0) += potential(state, state);
        }

        // export the energy mean
        Export::EigenMatrixDouble(std::string("POTENTIAL-ENERGY-MEAN_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the kinetic energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.kinetic_energy) {

        // create the matrix containing the energy
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, trajectory_data_vector.size());

        // fill the matrix with the energy
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) {
            int state = trajectory_data_vector.at(i).state(j); data_matrix.row(j)(i) = 0.5 * mass * trajectory_data_vector.at(i).velocity.row(j).squaredNorm();
        }

        // export the energy
        Export::EigenMatrixDouble(std::string("KINETIC-ENERGY_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix, time_variable);
    }

    // export the kinetic energy mean
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.kinetic_energy_mean) {

        // create the matrix containing the energy
        data_matrix = Eigen::MatrixXd::Zero(input.iterations + 1, 1);

        // fill the matrix with the energy
        for (int i = 0; i < input.trajectories; i++) for (int j = 0; j < input.iterations + 1; j++) {
            int state = trajectory_data_vector.at(i).state(j); data_matrix(j, 0) += 0.5 * mass * trajectory_data_vector.at(i).velocity.row(j).squaredNorm();
        }

        // export the energy mean
        Export::EigenMatrixDouble(std::string("KINETIC-ENERGY-MEAN_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix / input.trajectories, time_variable);
    }

    // export the hopping geometries
    if (input.data_export.hopping_geometry) {

        // create the data of the hopping geometries
        std::vector<Eigen::VectorXd> hopping_geometries;

        // fill the hopping geometry data
        for (int i = 0; i < input.trajectories; i++) for (size_t j = 0; j < trajectory_data_vector.at(i).hopping_geometry.size(); j++) {
            hopping_geometries.push_back(trajectory_data_vector.at(i).hopping_geometry.at(j));
        }

        // transform the hopping geometries to a matrix
        data_matrix = Eigen::MatrixXd::Zero(hopping_geometries.size(), hopping_geometries.front().rows()); for (size_t i = 0; i < hopping_geometries.size(); i++) data_matrix.row(i) = hopping_geometries.at(i).transpose();

        // export the energy mean
        Export::EigenMatrixDouble(std::string("HOPPING-GEOMETRIES_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix);
    }

    // export the hopping times
    if (input.data_export.hopping_time) {

        // create the data of the hopping times
        std::vector<double> hopping_times;

        // fill the hopping geometry data
        for (int i = 0; i < input.trajectories; i++) for (size_t j = 0; j < trajectory_data_vector.at(i).hopping_time.size(); j++) {
            hopping_times.push_back(trajectory_data_vector.at(i).hopping_time.at(j));
        }

        // transform the hopping geometries to a matrix
        data_matrix = Eigen::MatrixXd::Zero(hopping_times.size(), 1); for (size_t i = 0; i < hopping_times.size(); i++) data_matrix(i, 0) = hopping_times.at(i);

        // export the energy mean
        Export::EigenMatrixDouble(std::string("HOPPING-TIMES_") + algorithm + (input.adiabatic ? "-ADIABATIC" : "-DIABATIC") + ".mat", data_matrix);
    }
}

void Export::EigenMatrixDouble(const std::string& path, const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& independent_variable) {
    // create the matrix to export
    Eigen::MatrixXd export_matrix(matrix.rows(), matrix.cols() + independent_variable.cols());

    // copy the matrix and the independent variable to the export matrix
    if (independent_variable.size() > 0) export_matrix << independent_variable, matrix; else export_matrix = matrix;

    // open the output file and check for errors
    std::ofstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR WRITING");

    // write the dimensions to the header and set the formatting
    file << export_matrix.rows() << " " << export_matrix.cols() << "\n" << std::fixed << std::setprecision(16);

    // write the matrix by rows
    for (int i = 0; i < export_matrix.rows(); i++, file << "\n") for (int j = 0; j < export_matrix.cols(); j++) {
        file << std::setw(22) << export_matrix(i, j) << (j < export_matrix.cols() - 1 ? " " : "");
    }
}

void Export::TorchTensorDouble(const std::string& path, const torch::Tensor& tensor) {
    // open the input file and check for errors
    std::ofstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR WRITING");

    // write the dimensions to the header and set the formatting
    for (size_t i = 0; i < tensor.sizes().size(); i++) {file << tensor.sizes().at(i) << (i < tensor.sizes().size() - 1 ? " " : "");} file << "\n" << std::fixed << std::setprecision(16);

    // write 1st order tensor
    if (tensor.sizes().size() == 1) {
        for (int i = 0; i < tensor.sizes().at(0); i++) {
            file << std::setw(22) << tensor.index({i}).item().toDouble() << "\n";
        }
    }

    // write 2nd order tensor
    if (tensor.sizes().size() == 2) {
        for (int i = 0; i < tensor.sizes().at(0); i++) {
            for (int j = 0; j < tensor.sizes().at(1); j++) {
                file << std::setw(22) << tensor.index({i, j}).item().toDouble() << (j < tensor.sizes().at(1) - 1 ? " " : "");
            } file << "\n";
        }
    }

    // write 3rd order tensor
    if (tensor.sizes().size() == 3) {
        for (int i = 0; i < tensor.sizes().at(0); i++) {
            for (int j = 0; j < tensor.sizes().at(1); j++) {
                for (int k = 0; k < tensor.sizes().at(2); k++) {
                    file << std::setw(22) << tensor.index({i, j, k}).item().toDouble() << (j < tensor.sizes().at(1) - 1 ? " " : "");
                }
            } file << "\n";
        }
    }

    // write 4th order tensor
    if (tensor.sizes().size() == 4) {
        for (int i = 0; i < tensor.sizes().at(0); i++) {
            for (int j = 0; j < tensor.sizes().at(1); j++) {
                for (int k = 0; k < tensor.sizes().at(2); k++) {
                    for (int l = 0; l < tensor.sizes().at(3); l++) {
                        file << std::setw(22) << tensor.index({i, j, k, l}).item().toDouble() << (k < tensor.sizes().at(2) - 1 ? " " : "");
                    }
                } file << "\n";
            }
        }
    }

    // write 5th order tensor
    if (tensor.sizes().size() == 5) {
        for (int i = 0; i < tensor.sizes().at(0); i++) {
            for (int j = 0; j < tensor.sizes().at(1); j++) {
                for (int k = 0; k < tensor.sizes().at(2); k++) {
                    for (int l = 0; l < tensor.sizes().at(3); l++) {
                        for (int m = 0; m < tensor.sizes().at(4); m++) {
                            file << std::setw(22) << tensor.index({i, j, k, l, m}).item().toDouble() << (k < tensor.sizes().at(2) - 1 ? " " : "");
                        }
                    }
                } file << "\n";
            }
        }
    }
}

void Export::WavefunctionTrajectory(const Input::QuantumDynamics& input, const std::vector<QuantumDynamics::IterationData>& iteration_data, const Eigen::MatrixXd& grid, int state, bool imaginary) {
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
        Export::EigenMatrixDouble(std::string("PSI-DIABATIC_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, grid);
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
        Export::EigenMatrixDouble(std::string("PSI-ADIABATIC_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, grid);
    }

    // export the diabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.diabatic_density) {

        // create the container for the density matrix
        data_matrix = Eigen::MatrixXd(iteration_data.size(), input.potential.size() * input.potential.size());

        // fill the data matrix with diabatic density
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).density_diabatic.reshaped();

        // save the diabatic density matrix
        Export::EigenMatrixDouble(std::string("DENSITY-DIABATIC_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the adiabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.adiabatic_density) {

        // create the container for the density matrix
        data_matrix = Eigen::MatrixXd(iteration_data.size(), input.potential.size() * input.potential.size());

        // fill the data matrix with adiabatic density
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).density_adiabatic.reshaped();

        // save the adiabatic density matrix
        Export::EigenMatrixDouble(std::string("DENSITY-ADIABATIC_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the diabatic population matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.diabatic_population) {

        // create the container for the density matrix
        data_matrix = Eigen::MatrixXd(iteration_data.size(), input.potential.size());

        // fill the data matrix with diabatic populations
        for (size_t i = 0; i < iteration_data.size(); i++) for (size_t j = 0; j < input.potential.size(); j++) {
            data_matrix(i, j) = iteration_data.at(i).density_diabatic(j, j);
        }

        // save the diabatic population
        Export::EigenMatrixDouble(std::string("POPULATION-DIABATIC_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the adiabatic density matrix
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.adiabatic_population) {

        // create the container for the density matrix
        data_matrix = Eigen::MatrixXd(iteration_data.size(), input.potential.size());

        // fill the data matrix with adiabatic population matrix
        for (size_t i = 0; i < iteration_data.size(); i++) for (size_t j = 0; j < input.potential.size(); j++) {
            data_matrix(i, j) = iteration_data.at(i).density_adiabatic(j, j);
        }

        // save the adiabatic population matrix
        Export::EigenMatrixDouble(std::string("POPULATION-ADIABATIC_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the position
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.position) {

        // create the container for the position
        data_matrix = Eigen::MatrixXd(iteration_data.size(), iteration_data.front().position.size());

        // fill the data matrix with position of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).position;

        // save the wavefunction position
        Export::EigenMatrixDouble(std::string("POSITION_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the momentum
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.momentum) {

        // create the container for the momentum
        data_matrix = Eigen::MatrixXd(iteration_data.size(), iteration_data.front().momentum.size());

        // fill the data matrix with momentum of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix.row(i) = iteration_data.at(i).momentum;

        // save the wavefunction momentum
        Export::EigenMatrixDouble(std::string("MOMENTUM_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the total energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.total_energy) {

        // create the container for the energy
        data_matrix = Eigen::MatrixXd(iteration_data.size(), 1);

        // fill the data matrix with energy of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix(i) = iteration_data.at(i).total_energy;

        // save the wavefunction energy
        Export::EigenMatrixDouble(std::string("TOTAL-ENERGY_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the potential energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.potential_energy) {

        // create the container for the energy
        data_matrix = Eigen::MatrixXd(iteration_data.size(), 1);

        // fill the data matrix with energy of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix(i) = iteration_data.at(i).potential_energy;

        // save the wavefunction energy
        Export::EigenMatrixDouble(std::string("POTENTIAL-ENERGY_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the kinetic energy
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.kinetic_energy) {

        // create the container for the energy
        data_matrix = Eigen::MatrixXd(iteration_data.size(), 1);

        // fill the data matrix with energy of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix(i) = iteration_data.at(i).kinetic_energy;

        // save the wavefunction energy
        Export::EigenMatrixDouble(std::string("KINETIC-ENERGY_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }

    // export the autocorrelation function
    if (Eigen::VectorXd time_variable = Eigen::VectorXd::LinSpaced(input.iterations + 1, 0, input.time_step * input.iterations); input.data_export.acf) {

        // create the container for the autocorrelation function
        data_matrix = Eigen::MatrixXd(iteration_data.size(), 2);

        // fill the data matrix with autocorrelation function of the wavefunction
        for (size_t i = 0; i < iteration_data.size(); i++) data_matrix(i, 0) = iteration_data.at(i).acf.real(), data_matrix(i, 1) = iteration_data.at(i).acf.imag();

        // save the wavefunction autocorrelation function
        Export::EigenMatrixDouble(std::string("ACF_EXACT-") + (imaginary ? "IMAG" : "REAL") + "_"  + std::to_string(state + 1) + ".mat", data_matrix, time_variable);
    }
}

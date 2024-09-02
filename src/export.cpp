#include "export.h"

void Export::EigenMatrixDouble(const std::string& path, const Eigen::MatrixXd& matrix, const Eigen::VectorXd& independent_variable) {
    // create the matrix to export
    Eigen::MatrixXd export_matrix(matrix.rows(), matrix.cols() + (independent_variable.size() > 0 ? 1 : 0));

    // copy the matrix and the independent variable to the export matrix
    if (independent_variable.size() > 0) export_matrix << independent_variable, matrix; else export_matrix = matrix;

    // open the output file and write the dimensions to the header
    std::ofstream file(path); file << export_matrix.rows() << " " << export_matrix.cols() << "\n" << std::fixed << std::setprecision(16);

    // write the matrix by rows
    for (int i = 0; i < export_matrix.rows(); i++, file << "\n") for (int j = 0; j < export_matrix.cols(); j++) {
        file << std::setw(20) << export_matrix(i, j) << (j < export_matrix.cols() - 1 ? " " : "");
    }
}

void Export::WavefunctionTrajectory(const std::string& path, const std::vector<Wavefunction>& wavefunction_trajectory, const Eigen::VectorXd& independent_variable) {
    // create the container for the wavefunction trajectory
    Eigen::MatrixXd wavefunction_matrix(wavefunction_trajectory.front().get_data().rows(), 2 * wavefunction_trajectory.size() * wavefunction_trajectory.front().get_data().cols());

    // loop over all wavefunctions
    for (size_t i = 0; i < wavefunction_trajectory.size(); i++) {

        // extract the wavefunction data
        const Eigen::MatrixXcd& wavefunction_data = wavefunction_trajectory.at(i).get_data();

        // for each column of the data
        for (int j = 0; j < wavefunction_data.cols(); j++) {

            // extract the real and imaginary parts of the wavefunction column
            const Eigen::VectorXd& real_part = wavefunction_data.col(j).real();
            const Eigen::VectorXd& imag_part = wavefunction_data.col(j).imag();

            // copy the real and imaginary parts to the wavefunction matrix
            wavefunction_matrix.block(0, 2 * (i * wavefunction_data.cols() + j) + 0, wavefunction_data.rows(), 1) = real_part;
            wavefunction_matrix.block(0, 2 * (i * wavefunction_data.cols() + j) + 1, wavefunction_data.rows(), 1) = imag_part;
        }
    }

    // save the matrix
    Export::EigenMatrixDouble(path, wavefunction_matrix, independent_variable);
}

#include "export.h"

void Export::EigenMatrixDouble(const std::string& path, const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& independent_variable) {
    // create the matrix to export
    Eigen::MatrixXd export_matrix(matrix.rows(), matrix.cols() + independent_variable.cols());

    // copy the matrix and the independent variable to the export matrix
    if (independent_variable.size() > 0) export_matrix << independent_variable, matrix; else export_matrix = matrix;

    // open the output file and write the dimensions to the header
    std::ofstream file(path); file << export_matrix.rows() << " " << export_matrix.cols() << "\n" << std::fixed << std::setprecision(16);

    // write the matrix by rows
    for (int i = 0; i < export_matrix.rows(); i++, file << "\n") for (int j = 0; j < export_matrix.cols(); j++) {
        file << std::setw(20) << export_matrix(i, j) << (j < export_matrix.cols() - 1 ? " " : "");
    }
}

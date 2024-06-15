#include "linalg.h"
#include <fstream>
#include <iomanip>

Tensor<4> Eigen::Kron(const EigenMatrix<double>& A, const EigenTensor<double, 4>& B) {
    // define the tensor where the product will be stored
    EigenTensor<double, 4> C(B.dimension(0), B.dimension(1), A.rows() * B.dimension(2), A.cols() * B.dimension(3));

    // perform the Kronecker product
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            C.slice(array<Index, 4>({0, 0, i * B.dimension(2), j * B.dimension(3)}), array<Index, 4>({B.dimension(0), B.dimension(1), B.dimension(2), B.dimension(3)})) = A(i, j) * B;
        }
    }

    // return the Kronecker product
    return C;
}

Matrix Eigen::LoadMatrix(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and create the tensor
    int rows, cols; file >> rows >> cols; EigenMatrix<double> A(rows, cols);

    // read the tensor by dimensions, assign the values and return the matrix
    for (int i = 0; i < rows; i++) {for (int j = 0; j < cols; j++) file >> A(i, j);} return A;
}

Tensor<4> Eigen::LoadTensor(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and assign them to an array
    std::array<int, 4> dims; for (int i = 0; i < 4; i++) file >> dims.at(i);

    // create the tensor
    EigenTensor<double, 4> A(dims.at(0), dims.at(1), dims.at(2), dims.at(3));

    // read the tensor by dimensions, assign the values and return the tensor
    for (int i = 0; i < dims.at(0); i++) for (int j = 0; j < dims.at(1); j++) {
        for (int k = 0; k < dims.at(2); k++) for (int l = 0; l < dims.at(3); l++) file >> A(i, j, k, l);
    }

    // return the tensor
    return A;
}

void Eigen::Write(const std::string& path, const EigenTensor<double, 4>& A) {
    // open the output file and extract the dimensions
    std::ofstream file(path); auto dims = A.dimensions();

    // write the dimensions to the header and set the formatting
    for (size_t i = 0; i < dims.size(); i++) {file << dims.at(i) << (i < dims.size() - 1 ? " " : "");} file << "\n" << std::fixed << std::setprecision(14);

    // write the tensor by rows
    for (int i = 0; i < dims.at(0); i++) {
        for (int j = 0; j < dims.at(1); j++) {
            for (int k = 0; k < dims.at(2); k++) {
                for (int l = 0; l < dims.at(3); l++) {
                    file << std::setw(20) << A(i, j, k, l) << (k < dims.at(2) - 1 ? " " : "");
                }
            } file << "\n";
        }
    }
}

void Eigen::Write(const std::string& path, const EigenMatrix<double>& A) {
    // open the output file and write the dimensions to the header
    std::ofstream file(path); file << A.rows() << " " << A.cols() << "\n" << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < A.rows(); i++, file << "\n") for (int j = 0; j < A.cols(); j++) file << std::setw(20) << A(i, j) << (j < A.cols() - 1 ? " " : "");
}

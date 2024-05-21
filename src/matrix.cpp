#include "matrix.h"

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    Matrix<T> temp; temp.mat = mat + other.mat; return temp;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    Matrix<T> temp; temp.mat = mat + other.mat; return temp;
}

template <typename T>
Matrix<T> Matrix<T>::Load(const std::string& path) {
    // open the input file, read the dimensions and define the matrix
    std::ifstream file(path); int rows, cols; file >> cols >> rows; Matrix<T> mat(rows, cols);

    // read the matrix by rows
    for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) file >> mat(i, j);

    // return the matrix
    return Matrix<T>(mat);
}

template <typename T>
void Matrix<T>::save(const std::string& path) const {
    // open the output file
    std::ofstream file(path);

    // write the header and set the formatting
    file << mat.cols() << " " << mat.rows() << "\n" << std::fixed << std::setprecision(14);

    // write the matrix by rows
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            file << std::setw(24) << mat(i, j) << " ";
        } file << "\n";
    }
}

template class Matrix<double>;

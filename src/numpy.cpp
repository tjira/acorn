#include "numpy.h"

std::vector<std::vector<int>> Numpy::Combinations(int n, int k) {
    // get the number of combinations and create resulting array
    std::vector<std::vector<int>> combs(Factorial(n) / (Factorial(k) * Factorial(n - k)));

    // create the bitmask that will get permuted
    std::string bitmask(k, 1); bitmask.resize(n, 0);
 
    // generate the combinations
    for (size_t i = 0; i < combs.size(); i++) {
        for (int j = 0; j < n; j++) {
            if (bitmask[j]) combs.at(i).push_back(j);
        } std::prev_permutation(bitmask.begin(), bitmask.end());
    }

    // return the result
    return combs;
}

long Numpy::Factorial(long n) {
    return n > 1 ? n * Factorial(n - 1) : 1;
}

Matrix<> Numpy::Repeat(const Matrix<>& A, int count, int axis) {
    // create the result matrix
    Matrix<> B(axis ? A.rows() : count * A.rows(), axis ? count * A.cols() : A.cols());

    // throw error if axis is out of bounds
    if (axis != 0 && axis != 1) throw std::runtime_error("UNKNOWN AXIS IN MATRIX REPEAT");

    // repeat rows
    if (axis == 0) {
        for (int i = 0; i < A.rows(); i++) {
            for (int j = 0; j < count; j++) {
                B.row(i * count + j) = A.row(i);
            }
        }
    }

    // repeat columns
    else if (axis == 1) {
        for (int i = 0; i < A.cols(); i++) {
            for (int j = 0; j < count; j++) {
                B.col(i * count + j) = A.col(i);
            }
        }
    }

    // return result
    return B;
}

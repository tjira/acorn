#include "determinant.h"

std::vector<int> Common(const std::vector<int>& first, const std::vector<int>& second) {
    // create the container
    std::vector<int> common;

    // fill the container
    for (size_t i = 0; i < first.size(); i++) if (first.at(i) == second.at(i)) common.push_back(first.at(i));

    // return the container
    return common;
}

bool VectorContains(const std::vector<int>& v, const int& e) {return std::find(v.begin(), v.end(), e) != v.end();}

std::vector<int> Unique(const std::vector<int>& first, const std::vector<int>& second) {
    // create the container
    std::vector<int> unique;

    // fill the container
    for (size_t i = 0; i < second.size(); i++) if (!VectorContains(first, second.at(i))) unique.push_back(second.at(i));
    for (size_t i = 0; i < first.size(); i++) if (!VectorContains(second, first.at(i))) unique.push_back(first.at(i));

    // return the container
    return unique;
}

Determinant::Determinant(int norb, const std::vector<int>& a, const std::vector<int>& b) : a(a), b(b), norb(norb) {}

Determinant::Determinant(int norb, int nocca, int noccb) : a(nocca), b(noccb), norb(norb) {
    std::iota(a.begin(), a.end(), 0), std::iota(b.begin(), b.end(), 0);
}

std::tuple<Determinant, int> Determinant::align(const Determinant& second) const {
    // define the temporary determinant and swaps variable
    Determinant first(*this); int swaps = 0;
    
    // align the alpha electrons
    for (size_t i = 0; i < first.a.size(); i++) {
        if (first.a.at(i) != second.a.at(i)) {
            for (size_t j = 0; j < a.size(); j++) {
                if (first.a.at(i) == second.a.at(j)) {
                    std::swap(first.a.at(i), first.a.at(j)), swaps++;
                }
            }
        }
    }

    // align the beta electrons
    for (size_t i = 0; i < first.b.size(); i++) {
        if (first.b.at(i) != second.b.at(i)) {
            for (size_t j = 0; j < b.size(); j++) {
                if (first.b.at(i) == second.b.at(j)) {
                    std::swap(first.b.at(i), first.b.at(j)), swaps++;
                }
            }
        }
    }

    // return the aligned determinant
    return std::tuple<Determinant, int>{first, swaps};
}

int Determinant::differences(const Determinant& second) const {
    // define the count
    int count = 0;

    // count the differences
    for (size_t i = 0; i < a.size(); i++) if (a.at(i) != second.a.at(i)) count++;
    for (size_t i = 0; i < b.size(); i++) if (b.at(i) != second.b.at(i)) count++;

    // return the count
    return count;
}

std::vector<Determinant> Determinant::full() const {
    // define the vector of determinants
    std::vector<Determinant> full;

    // generate all the combinations
    for (const std::vector<int>& alpha : Numpy::Combinations(norb, a.size())) {
        for (const std::vector<int>& beta : Numpy::Combinations(norb, b.size())) {
            full.push_back(Determinant(norb, alpha, beta));
        }
    }

    // return vector
    return full;
}

#include <iostream>

double Determinant::hamilton(const Determinant& second, const Matrix<>& Hms, const Tensor<4>& Jms) const {
    // align this determinant, number of swaps and differences and the matrix element
    auto[first, swaps] = align(second); int diff = first.differences(second); double elem = 0;
    std::vector<int> firstso = first.spinorbitals(), secondso = second.spinorbitals();

    // get the common and unique  spinorbitals
    std::vector<int> common = Common(firstso, secondso), unique = Unique(firstso, secondso);

    // assign the matrix element
    if (diff == 0) {
        for (int so : firstso) elem += Hms(so, so);
        for (size_t i = 0; i < firstso.size() - 1; i++) {
            for (size_t j = i + 1; j < firstso.size(); j++) {
                elem += Jms(firstso.at(i), firstso.at(j), firstso.at(i), firstso.at(j));
            }
        }
    } else if (diff == 1) {
        elem += Hms(unique.at(0), unique.at(1));
        for (int so : common) {
            elem += Jms(unique.at(0), so, unique.at(1), so);
        }
    } else if (diff == 2) {
        elem = Jms(unique.at(0), unique.at(3), unique.at(2), unique.at(1));
    }

    // return element
    return std::pow(-1, swaps) * elem;
}

std::vector<int> Determinant::spinorbitals() const {
    // create the vector to hold all spinorbitals
    std::vector<int> spinorbitals(a.size() + b.size());

    // fill the spinorbitals
    for (size_t i = 0; i < a.size(); i++) spinorbitals.at(i) = 2 * a.at(i);
    for (size_t i = 0; i < b.size(); i++) spinorbitals.at(a.size() + i) = 2 * b.at(i) + 1;

    // return the spinorbital vector
    return spinorbitals;
}

#include "determinant.h"

bool VectorContains(const std::vector<int>& v, const int& e) {return std::find(v.begin(), v.end(), e) != v.end();}

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

double Determinant::hamilton(const Determinant& second, const Matrix<>& Hms, const Tensor<4>& Jms) const {
    // align this determinant and define differences and element variable
    auto[first, swaps] = align(second); int diff = 0; double elem = 0;

    // define all needed vectors of spinorbitals
    std::vector<int> firstso(a.size() + b.size()), secondso(a.size() + b.size()), common, unique;

    // calculate the number of differences
    for (size_t i = 0; i < a.size(); i++) if (first.a.at(i) != second.a.at(i)) diff++;
    for (size_t i = 0; i < b.size(); i++) if (first.b.at(i) != second.b.at(i)) diff++;

    // fill the spinorbitals of first and second determinant
    for (size_t i = 0; i < a.size(); i++) secondso.at(i) = 2 * second.a.at(i), secondso.at(second.a.size() + i) = 2 * second.b.at(i) + 1;
    for (size_t i = 0; i < a.size(); i++) firstso.at(i) = 2 * first.a.at(i), firstso.at(first.a.size() + i) = 2 * first.b.at(i) + 1;

    // fill the common spinorbitals vector
    for (size_t i = 0; i < firstso.size(); i++) if (firstso.at(i) == secondso.at(i)) common.push_back(firstso.at(i));

    // fill the unique spinorbitals vector
    for (size_t i = 0; i < firstso.size(); i++) if (!VectorContains(firstso, secondso.at(i))) unique.push_back(secondso.at(i));
    for (size_t i = 0; i < firstso.size(); i++) if (!VectorContains(secondso, firstso.at(i))) unique.push_back(firstso.at(i));

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

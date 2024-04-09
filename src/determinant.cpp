#include "determinant.h"

Determinant::Determinant(int norb, const std::vector<int>& a, const std::vector<int>& b) : a(a), b(b), norb(norb) {}

Determinant::Determinant(int norb, int nocca, int noccb) : a(nocca), b(noccb), norb(norb) {
    std::iota(a.begin(), a.end(), 0), std::iota(b.begin(), b.end(), 0);
}

std::tuple<Determinant, int> Determinant::align(const Determinant& det2) const {
    // define the temporary determinant and swaps variable
    Determinant det1(*this); int swaps = 0;
    
    // align the alpha electrons
    for (size_t i = 0; i < det1.a.size(); i++) {
        if (det1.a.at(i) != det2.a.at(i)) {
            for (size_t j = 0; j < a.size(); j++) {
                if (det1.a.at(i) == det2.a.at(j)) {
                    std::swap(det1.a.at(i), det1.a.at(j)), swaps++;
                }
            }
        }
    }

    // align the beta electrons
    for (size_t i = 0; i < det1.b.size(); i++) {
        if (det1.b.at(i) != det2.b.at(i)) {
            for (size_t j = 0; j < b.size(); j++) {
                if (det1.b.at(i) == det2.b.at(j)) {
                    std::swap(det1.b.at(i), det1.b.at(j)), swaps++;
                }
            }
        }
    }

    // return the aligned determinant
    return std::tuple<Determinant, int>{det1, swaps};
}

std::vector<Determinant> Determinant::excitations(const std::vector<int>& excs) const {
    // check if the excitations are valid
    if (excs.size() > 2 || std::count_if(excs.begin(), excs.end(), [](int exc) {return exc > 2;}) > 0) {
        throw std::runtime_error("PROVIDED EXCITATIONS NOT IMPLEMENTED IN CI");
    }

    // define the vector of excitations
    std::vector<Determinant> excitations = {*this};

    // single excitations
    if (std::find(excs.begin(), excs.end(), 1) != excs.end()) {
        for (int i = 0; i < a.size(); i++) {
            for (int j = a.size(); j < norb; j++) {
                std::vector<int> det = this->a; det.at(i) = j; excitations.push_back(Determinant(norb, det, b));
            }
        }
        for (int i = 0; i < b.size(); i++) {
            for (int j = b.size(); j < norb; j++) {
                std::vector<int> det = this->b; det.at(i) = j; excitations.push_back(Determinant(norb, a, det));
            }
        }
    }

    // double excitations
    if (std::find(excs.begin(), excs.end(), 2) != excs.end()) {
        for (int i = 0; i < a.size(); i++) {
            for (int j = a.size(); j < norb; j++) {
                for (int k = 0; k < b.size(); k++) {
                    for (int l = b.size(); l < norb; l++) {
                        std::vector<int> deta = this->a, detb = this->b; deta.at(i) = j, detb.at(k) = l; excitations.push_back(Determinant(norb, deta, detb));
                    }
                }
            }
        }
        for (int i = 0; i < a.size(); i++) {
            for (int j = a.size(); j < norb; j++) {
                for (int k = i + 1; k < a.size(); k++) {
                    for (int l = j + 1; l < norb; l++) {
                        std::vector<int> det = this->a; det.at(i) = j, det.at(k) = l; excitations.push_back(Determinant(norb, det, b));
                    }
                }
            }
        }
        for (int i = 0; i < b.size(); i++) {
            for (int j = b.size(); j < norb; j++) {
                for (int k = i + 1; k < b.size(); k++) {
                    for (int l = j + 1; l < norb; l++) {
                        std::vector<int> det = this->b; det.at(i) = j, det.at(k) = l; excitations.push_back(Determinant(norb, a, det));
                    }
                }
            }
        }
    }

    // return the excitations
    return excitations;
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

double Determinant::hamilton(const Determinant& det2, const Matrix<>& Hms, const Tensor<4>& Jms) const {
    // align this determinant and define differences and element variable
    auto[det1, swaps] = align(det2); int diff = 0; double elem = 0;

    // define all needed vectors of spinorbitals
    std::vector<int> so1(a.size() + b.size()), so2(a.size() + b.size()), common, unique;

    // calculate the number of differences
    for (size_t i = 0; i < a.size(); i++) if (det1.a.at(i) != det2.a.at(i)) diff++;
    for (size_t i = 0; i < b.size(); i++) if (det1.b.at(i) != det2.b.at(i)) diff++;

    // fill the spinorbitals of det1 and det2 determinant
    for (size_t i = 0; i < a.size(); i++) so2.at(i) = 2 * det2.a.at(i), so2.at(det2.a.size() + i) = 2 * det2.b.at(i) + 1;
    for (size_t i = 0; i < a.size(); i++) so1.at(i) = 2 * det1.a.at(i), so1.at(det1.a.size() + i) = 2 * det1.b.at(i) + 1;

    // fill the common spinorbitals vector
    for (size_t i = 0; i < so1.size(); i++) if (so1.at(i) == so2.at(i)) common.push_back(so1.at(i));

    // fill the unique spinorbitals vector
    for (size_t i = 0; i < so1.size(); i++) if (!VectorContains(so1, so2.at(i))) unique.push_back(so2.at(i));
    for (size_t i = 0; i < so1.size(); i++) if (!VectorContains(so2, so1.at(i))) unique.push_back(so1.at(i));

    // assign the matrix element according to Slater-Condon rules
    if (diff == 0) {
        for (int so : so1) elem += Hms(so, so);
        for (size_t i = 0; i < so1.size(); i++) {
            for (size_t j = 0; j < so1.size(); j++) {
                elem += 0.5 * (Jms(so1.at(i), so1.at(i), so1.at(j), so1.at(j)) - Jms(so1.at(i), so1.at(j), so1.at(j), so1.at(i)));
            }
        }
    } else if (diff == 1) {
        elem += Hms(unique.at(0), unique.at(1));
        for (int so : common) {
            elem += Jms(unique.at(0), unique.at(1), so, so) - Jms(unique.at(0), so, so, unique.at(1));
        }
    } else if (diff == 2) {
        elem = Jms(unique.at(0), unique.at(2), unique.at(1), unique.at(3)) - Jms(unique.at(0), unique.at(3), unique.at(1), unique.at(2));
    }

    // return element
    return std::pow(-1, swaps) * elem;
}

#pragma once

#include "numpy.h"

class Determinant {
public:
    // constructors
    Determinant(int norb, int nocca, int noccb); Determinant(int norb, const std::vector<int>& a, const std::vector<int>& b);

    // excitation generator and aligner
    std::tuple<Determinant, int> align(const Determinant& second) const; std::vector<Determinant> excitations(const std::vector<int>& excs) const; std::vector<Determinant> full() const;

    // integrals and functions
    double hamilton(const Determinant& second, const Matrix<>& Hms, const Tensor<4>& Jms) const;

private:
    std::vector<int> a, b; int norb;
};

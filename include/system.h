#pragma once

#include "matrix.h"

class System {
public:
    System(const std::string& path);

    // property calculators
    double nuclearRepulsion() const;

    // private variable getters
    Matrix<> atoms() const; Matrix<> coords() const; int nocc() const;

private:
    Matrix<> AN, R;
};

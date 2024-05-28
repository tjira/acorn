#pragma once

#include "linalg.h"

class System {
public:
    System(const std::string& path);

    // property calculators
    double nuclearRepulsion() const; int nocc() const;

private:
    EigenMatrix<> AN, R;
};

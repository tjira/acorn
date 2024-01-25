#pragma once

#include "expression.h"

class ModelSystem {
    friend class ModelSolver;
public:
    // constructors and destructors
    ModelSystem(int m, const std::vector<std::vector<std::string>>& potential, const std::vector<std::vector<double>> limits, int ngrid);

    // static methods
    void static SaveWavefunction(const std::string& fname, const Vector<>& r, const std::vector<Vector<std::complex<double>>>& wfns, const std::vector<double>& energy);

private:
    const std::vector<std::vector<double>> limits; const int ngrid, m;
    const std::vector<std::vector<std::string>> potential;
};

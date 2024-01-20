#pragma once

#include "expression.h"

class ModelSystem {
    friend class ModelSolver;
public:
    // constructors and destructors
    ModelSystem(const std::string& potential, const std::vector<std::vector<double>> limits, int ngrid);

    // static methods
    void static SaveWavefunction(const std::string& fname, const Vector<>& r, const std::vector<std::vector<Vector<std::complex<double>>>>& wfns, const Vector<>& energy);

private:
    const std::vector<std::vector<double>> limits;
    const std::string potential; const int ngrid;
};

#pragma once

#include "expression.h"

class ModelSystem {
    friend class ModelSolver;
public:
    // constructors and destructors
    ModelSystem(int m, const std::vector<std::vector<std::string>>& potential, const std::vector<std::string>& variables, const std::vector<double> limits, int ngrid);

    // static methods
    void static SaveWavefunction(const std::filesystem::path& fname, const Vector<>& x, const Vector<>& y, const std::vector<Matrix<std::complex<double>>>& wfns, const std::vector<double>& energy);
    void static SaveWavefunction(const std::filesystem::path& fname, const Vector<>& x, const std::vector<Matrix<std::complex<double>>>& wfns, const std::vector<double>& energy);

    // getters
    int mass() const {return m;} std::vector<std::string> vars() const {return variables;}

private:
    const std::vector<std::vector<std::string>> potential;
    const std::vector<double> limits; const int ngrid, m;
    const std::vector<std::string> variables;
};

#pragma once

#include "default.h"
#include "eigen.h"
#include "timer.h"

inline std::unordered_map<int, std::string> an2sm = {
    { 1,  "H"},
    { 6,  "C"},
    { 7,  "N"},
    { 8,  "O"},
    { 9,  "F"},
    {17, "Cl"}
};

class System {
public:
    System(std::ifstream& stream, std::string basis, int charge = 0, int multi = 1);
    double repulsion() const; libint2::BasisSet getShells() const {return shells;}
    void save(std::string fname, std::ios::openmode mode = std::ios::out) const;
    int nocc() const {return electrons / 2;} void move(const Matrix<>& dir);
    std::vector<libint2::Atom> getAtoms() const {return atoms;}

private:
    std::vector<libint2::Atom> atoms;
    int electrons, charge, multi;
    libint2::BasisSet shells;
    std::string basis;
};

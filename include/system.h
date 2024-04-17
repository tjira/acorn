#pragma once

#include "eigen.h"
#include "table.h"

namespace libint2 {
    class BasisSet; int major(); int minor(); int micro();
} struct Atom {int atomic_number; double x, y, z;};

class System {
public:
    // constructors and destructors
    ~System(); System(const System& system); System(std::ifstream& stream, std::string basis, int charge = 0, int multi = 1);

    // printer operator and electron getter
    friend std::ostream& operator<<(std::ostream& os, const System& system); int getElectrons() const {return electrons;}

    // atom container getter and charge, multi and basis getters
    int getCharge() const {return charge;} int getMulti() const {return multi;} std::string getBasis() const {return lnxbasis;}
    template <typename T = Atom> std::vector<T> getAtoms() const {return *reinterpret_cast<const std::vector<T>*>(&atoms);}

    // shell container and getter, system mover, number of occupied orbitals getter
    libint2::BasisSet getShells() const; void move(const Matrix<>& dir); int nocc() const {return electrons / 2;}

    // molecule exporter and nuclear repulsion getter
    void save(std::filesystem::path fname, std::ios::openmode mode = std::ios::out) const; double repulsion() const; Matrix<> drepulsion() const;

private:
    int electrons, charge, multi; std::string basis, lnxbasis;
    libint2::BasisSet* shells; std::vector<Atom> atoms;
};

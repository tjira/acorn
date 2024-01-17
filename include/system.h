#pragma once

#include "eigen.h"
#include "table.h"

namespace libint2 {
    class BasisSet; int major(); int minor(); int micro();
} struct Atom {int atomic_number; double x, y, z;};

class System {
public:
    // constructors and destructors
    System(std::ifstream& stream, std::string basis, int charge = 0, int multi = 1);

    // atom container getter
    template <typename T = Atom> std::vector<T> getAtoms() const {return *reinterpret_cast<const std::vector<T>*>(&atoms);}

    // shell container and getter, molecule mover, and number of occupied orbitals getter
    libint2::BasisSet getShells() const; void move(const Matrix<>& dir); int nocc() const {return electrons / 2;}

    // molecule exporter and nuclear repulsion getter
    void save(std::string fname, std::ios::openmode mode = std::ios::out) const; double repulsion() const;

private:
    int electrons, charge, multi;
    std::vector<Atom> atoms;
    std::string basis;
};

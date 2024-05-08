#include "population.h"

// include libint
#include "libint.h"

Vector<> Population::Mulliken(const System& system, const Integrals& ints, const Matrix<>& D) {
    // define the resulting vector and the DS matrix
    Matrix<> q = Vector<>::Zero(system.getAtoms().size()), DS = D.cwiseProduct(ints.S);

    // add the contribution from molecular orbitals
    for (size_t i = 0, j = 0; i < system.getShells().size(); i++) {
        for (size_t k = 0; k < system.getShells().at(i).size(); k++) {
            q(system.getShells().shell2atom(system.getAtoms<libint2::Atom>()).at(i)) -= DS.colwise().sum()(j + k);
        } j += system.getShells().at(i).size();
    }

    // add the contribution from the atomic number
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        q(i) += system.getAtoms().at(i).atomic_number;
    }

    // return the population vector
    return q;
}

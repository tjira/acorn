#include "system.h"

System::System(std::ifstream& stream, std::string basis, int charge, int multi) : electrons(0), charge(charge), multi(multi), basis(basis) {
    // check for the input geometry file existence
    if (!stream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // replace the basis placeholders
    std::replace(this->basis.begin(), this->basis.end(), '*', 's'), std::replace(this->basis.begin(), this->basis.end(), '+', 'p');

    // read the input geometry
    atoms = libint2::read_dotxyz(stream);

    // assign shells and calculate the number of electrons
    shells = libint2::BasisSet(basis, atoms, true); electrons -= charge;
    for (const auto& atom : atoms) electrons += atom.atomic_number;
}

void System::move(const Matrix<>& dir) {
    // shift the atoms
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms.at(i).x += dir(i, 0);
        atoms.at(i).y += dir(i, 1);
        atoms.at(i).z += dir(i, 2);
    }

    // shift the coords and recreate the basis
    shells = libint2::BasisSet(basis, atoms, true);
}

double System::repulsion() const {
    // create nuclear repulsion variable
    double repulsion = 0;

    for (size_t i = 0; i < atoms.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            // extract the atoms at position i and j
            libint2::Atom a = atoms.at(i), b = atoms.at(j);

            // calculate the vector
            double x = a.x - b.x, y = a.y - b.y, z = a.z - b.z;

            // add the value to the repulsion
            repulsion += a.atomic_number * b.atomic_number / std::sqrt(x * x + y * y + z * z);
        }
    }

    // return the repulsion
    return repulsion;
}

void System::save(std::string fname, std::ios::openmode mode) const {
    // open the file, write number of atoms and set output stream flags
    std::ofstream file(fname, mode); file << atoms.size() << "\n" << fname;
    file << std::fixed << std::setprecision(14) << "\n";

    // print all the atom coordinates
    for (size_t i = 0; i < atoms.size(); i++) {
        file << an2sm.at(atoms.at(i).atomic_number) << " "
             << std::setw(20) << BOHR2A * atoms.at(i).x << " "
             << std::setw(20) << BOHR2A * atoms.at(i).y << " "
             << std::setw(20) << BOHR2A * atoms.at(i).z << std::endl;
    }
}

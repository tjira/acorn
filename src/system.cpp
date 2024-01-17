#include "system.h"

// include the libint
#include <libint2.hpp>

// define the functions returning libint version
int libint2::major() {return LIBINT_MAJOR_VERSION;}
int libint2::minor() {return LIBINT_MINOR_VERSION;}
int libint2::micro() {return LIBINT_MICRO_VERSION;}

// define the shell getter
libint2::BasisSet System::getShells() const {return libint2::BasisSet(basis, *reinterpret_cast<const std::vector<libint2::Atom>*>(&atoms), true);}

System::System(std::ifstream& stream, std::string basis, int charge, int multi) : electrons(0), charge(charge), multi(multi), basis(basis) {
    // check for the input geometry file existence
    if (!stream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // replace the basis placeholders
    std::replace(this->basis.begin(), this->basis.end(), '*', 's'), std::replace(this->basis.begin(), this->basis.end(), '+', 'p');

    // read the input geometry
    std::vector<libint2::Atom> libintatoms = libint2::read_dotxyz(stream);
    atoms = *reinterpret_cast<const std::vector<Atom>*>(&libintatoms);

    // calculate the number of electrons
    electrons -= charge; for (const auto& atom : atoms) electrons += atom.atomic_number;
}


void System::move(const Matrix<>& dir) {
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms.at(i).x += dir(i, 0), atoms.at(i).y += dir(i, 1), atoms.at(i).z += dir(i, 2);
    }
}

double System::repulsion() const {
    // create the variable
    double repulsion = 0;

    for (size_t i = 0; i < atoms.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            // extract the atoms at position i and j and calculate their distance
            auto a = atoms.at(i), b = atoms.at(j); double x = a.x - b.x, y = a.y - b.y, z = a.z - b.z;

            // add the value to the repulsion
            repulsion += a.atomic_number * b.atomic_number / std::sqrt(x * x + y * y + z * z);
        }
    }

    // return the repulsion
    return repulsion;
}

void System::save(std::string fname, std::ios::openmode mode) const {
    // open the file, write number of atoms and set output stream flags
    std::ofstream file(fname, mode); file << atoms.size() << "\n";
    file << fname << std::fixed << std::setprecision(14) << "\n";

    // print all the atom coordinates
    for (size_t i = 0; i < atoms.size(); i++) {
        file << an2sm.at(atoms.at(i).atomic_number) << " ";
        file << std::setw(20) << BOHR2A * atoms.at(i).x << " ";
        file << std::setw(20) << BOHR2A * atoms.at(i).y << " ";
        file << std::setw(20) << BOHR2A * atoms.at(i).z << "\n";
    }
}

template std::vector<libint2::Atom> System::getAtoms() const;
template std::vector<Atom> System::getAtoms() const;

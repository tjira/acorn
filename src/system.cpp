#include "system.h"

// include the libint
#include <libint2.hpp>

// define the functions returning libint version
int libint2::major() {return LIBINT_MAJOR_VERSION;}
int libint2::minor() {return LIBINT_MINOR_VERSION;}
int libint2::micro() {return LIBINT_MICRO_VERSION;}

// define the shell getter and the destructor
libint2::BasisSet System::getShells() const {return *shells;} System::~System() {delete shells;}

System::System(const System& system) : electrons(system.electrons), charge(system.charge), multi(system.multi), basis(system.basis), lnxbasis(system.lnxbasis), atoms(system.atoms) {
    shells = new libint2::BasisSet(basis, *reinterpret_cast<const std::vector<libint2::Atom>*>(&atoms));
}

System::System(std::ifstream& stream, std::string basis, int charge, int multi) : electrons(0), charge(charge), multi(multi), basis(basis), lnxbasis(basis) {
    // check for the input geometry file existence
    if (!stream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // throw an error if impossible combination of charge and multiplicity
    if (std::abs(charge) % 2 == 0 && multi % 2 == 0) {
        throw std::runtime_error("MOLECULE CAN'T HAVE AN EVEN CHARGE AND MULTIPLICITY AT THE SAME TIME.");
    } else if (std::abs(charge) % 2 == 1 && multi % 2 == 1) {
        throw std::runtime_error("MOLECULE CAN'T HAVE AN ODD CHARGE AND MULTIPLICITY AT THE SAME TIME.");
    }

    // replace the basis placeholders
    std::replace(this->basis.begin(), this->basis.end(), '*', 's'), std::replace(this->basis.begin(), this->basis.end(), '+', 'p');

    // read the input geometry
    std::vector<libint2::Atom> libintatoms = libint2::read_dotxyz(stream);
    atoms = *reinterpret_cast<const std::vector<Atom>*>(&libintatoms);

    // create the shell container
    shells = new libint2::BasisSet(basis, *reinterpret_cast<const std::vector<libint2::Atom>*>(&atoms));

    // calculate the number of electrons
    electrons -= charge; for (const auto& atom : atoms) electrons += atom.atomic_number;
}

std::ostream& operator<<(std::ostream& os, const System& system) {
    // print the number of atoms and the basis
    os << system.atoms.size() << "\n" << system.lnxbasis << "\n";

    // print all the atom coordinates
    for (size_t i = 0; i < system.atoms.size(); i++) {
        os << std::setw(2) << an2sm.at(system.atoms.at(i).atomic_number) << " ";
        os << std::setw(20) << BOHR2A * system.atoms.at(i).x << " ";
        os << std::setw(20) << BOHR2A * system.atoms.at(i).y << " ";
        os << std::setw(20) << BOHR2A * system.atoms.at(i).z;
        if (i < system.atoms.size() - 1) os << "\n";
    }

    // return the output stream
    return os;
}

void System::move(const Matrix<>& dir) {
    // move the system
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms.at(i).x += dir(i, 0), atoms.at(i).y += dir(i, 1), atoms.at(i).z += dir(i, 2);
    }

    // create the shell container
    delete shells; shells = new libint2::BasisSet(basis, *reinterpret_cast<const std::vector<libint2::Atom>*>(&atoms));
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

    // check for correct file output
    if (!file.good()) throw std::runtime_error("ERROR WHEN OPENING THE OUTPUT FILE");

    // print all the atom coordinates
    for (size_t i = 0; i < atoms.size(); i++) {
        file << std::setw(2) << an2sm.at(atoms.at(i).atomic_number) << " ";
        file << std::setw(20) << BOHR2A * atoms.at(i).x << " ";
        file << std::setw(20) << BOHR2A * atoms.at(i).y << " ";
        file << std::setw(20) << BOHR2A * atoms.at(i).z << "\n";
    }
}

template std::vector<libint2::Atom> System::getAtoms() const;
template std::vector<Atom> System::getAtoms() const;

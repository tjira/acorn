#include "system.h"

System::System(const Input::System& input) : charge(input.charge), multi(input.multiplicity), basis(input.basis) {
    // open the file stream, define the periodic table and initialize container variables
    std::ifstream system_file_stream(input.path); std::vector<std::string> ptable; std::string line, sm; int natoms;

    // throw an error if impossible combination of charge and multiplicity
    if (std::abs(charge) % 2 == 0 && multi % 2 == 0) throw std::runtime_error("MOLECULE CAN'T HAVE AN EVEN CHARGE AND MULTIPLICITY AT THE SAME TIME.");
    if (std::abs(charge) % 2 == 1 && multi % 2 == 1) throw std::runtime_error("MOLECULE CAN'T HAVE AN ODD CHARGE AND MULTIPLICITY AT THE SAME TIME." );

    // check if the file stream is good
    if (!system_file_stream.good()) throw std::runtime_error("COULD NOT OPEN THE `" + input.path + "` FILE");

    // extract the number of atoms, skip the first two lines and initialize the atomic number and position matrices
    std::getline(system_file_stream, line); std::stringstream(line) >> natoms; std::getline(system_file_stream, line); atomic_numbers.resize(natoms); coordinates = torch::zeros({natoms, 3}, torch::kDouble);

    // extract the atomic numbers and positions
    for (int i = 0; i < natoms; i++) {

        // extract the atomic number and coordinates
        double x, y, z; std::getline(system_file_stream, line), std::stringstream(line) >> sm >> x >> y >> z;

        // assign the atomic number and position
        coordinates.index_put_({i, 0}, x), coordinates.index_put_({i, 1}, y), coordinates.index_put_({i, 2}, z), atomic_numbers.at(i) = sm2an.at(sm);
    }
}

int System::basis_functions() const {
    return get_shells().nbf();
}

int System::electrons() const {
    return std::accumulate(atomic_numbers.begin(), atomic_numbers.end(), 0) - charge;
}

std::vector<libint2::Atom> System::get_atoms() const {
    // create the vector of atoms
    std::vector<libint2::Atom> atoms;

    // loop over the atomic numbers and coordinates
    for (int i = 0; i < (int)atomic_numbers.size(); i++) {

        // extract the coordinates
        double x = coordinates.index({i, 0}).item<double>() * ANGSTROM_TO_BOHR;
        double y = coordinates.index({i, 1}).item<double>() * ANGSTROM_TO_BOHR;
        double z = coordinates.index({i, 2}).item<double>() * ANGSTROM_TO_BOHR;

        // add the atom to the vector
        atoms.push_back({atomic_numbers.at(i), x, y, z});
    }

    // return the atoms
    return atoms;
}

std::string System::get_basis() const {
    std::string uppercase_basis = basis; std::transform(uppercase_basis.begin(), uppercase_basis.end(), uppercase_basis.begin(), ::toupper); return uppercase_basis;
}

int System::get_multi() const {
    return multi;
}

libint2::BasisSet System::get_shells() const {
    // define the basis file
    std::string basis_file = basis;

    // replace all the * with s and + with p
    std::replace(basis_file.begin(), basis_file.end(), '*', 's');
    std::replace(basis_file.begin(), basis_file.end(), '+', 'p');

    // replace the augmented keyword
    basis_file = std::regex_replace(basis_file, std::regex("aug-"), "augmented-");

    // return the basis set
    return libint2::BasisSet(basis_file, get_atoms());
}

double System::nuclear_repulsion() const {
    // create the variable
    double nuclear_repulsion_value = 0;

    // calculate the repulsion
    for (int i = 0; i < (int)atomic_numbers.size(); i++) for (int j = 0; j < i; j++) {
        double d = (coordinates.index({i, None}) - coordinates.index({j, None})).norm().item<double>(); nuclear_repulsion_value += atomic_numbers.at(i) * atomic_numbers.at(j) / d / ANGSTROM_TO_BOHR;
    }

    // return the nuclear-nuclear repulsion of the system
    return nuclear_repulsion_value;
}

torch::Tensor System::nuclear_repulsion_d1() const {
    // create the variable
    torch::Tensor nuclear_repulsion_value = torch::zeros_like(coordinates);

    // calculate the repulsion derivative
    for (int i = 0; i < (int)atomic_numbers.size(); i++) for (int j = 0; j < (int)atomic_numbers.size(); j++) {

        // skip the same atom
        if (i == j) continue;

        // extract the vector
        torch::Tensor vector = (coordinates.index({i, None}) - coordinates.index({j, None})).squeeze();

        // calculate the derivative
        nuclear_repulsion_value.index({i, 0}) -= vector.index({0}).item<double>() * atomic_numbers.at(i) * atomic_numbers.at(j) / std::pow(vector.norm().item<double>(), 3) / std::pow(ANGSTROM_TO_BOHR, 2);
        nuclear_repulsion_value.index({i, 1}) -= vector.index({1}).item<double>() * atomic_numbers.at(i) * atomic_numbers.at(j) / std::pow(vector.norm().item<double>(), 3) / std::pow(ANGSTROM_TO_BOHR, 2);
        nuclear_repulsion_value.index({i, 2}) -= vector.index({2}).item<double>() * atomic_numbers.at(i) * atomic_numbers.at(j) / std::pow(vector.norm().item<double>(), 3) / std::pow(ANGSTROM_TO_BOHR, 2);
    }

    // return the nuclear-nuclear repulsion of the system
    return nuclear_repulsion_value;
}

int System::virtual_spinorbitals() const {
    return 2 * basis_functions() - electrons();
}

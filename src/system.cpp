#include "system.h"

System::System(const Input::System& input) : basis(input.basis) {
    // open the file stream, define the periodic table and initialize container variables
    std::ifstream system_file_stream(input.path); std::vector<std::string> ptable; std::string line, sm; int natoms;

    // check if the file stream is good
    if (!system_file_stream.good()) throw std::runtime_error("COULD NOT OPEN THE `" + input.path + "` FILE");

    // fill the periodic table
    for (std::stringstream periodic_table_stream(PERIODIC_TABLE); periodic_table_stream.good(); ptable.push_back(""), periodic_table_stream >> ptable.back()) {}

    // define the symbol to atomic number map
    auto symbol_to_atomic_number_map = [&ptable](const std::string& sm) {return std::find(ptable.begin(), ptable.end(), sm) - ptable.begin();};

    // extract the number of atoms, skip the first two lines and initialize the atomic number and position matrices
    std::getline(system_file_stream, line); std::stringstream(line) >> natoms; std::getline(system_file_stream, line); atomic_numbers.resize(natoms); coordinates = torch::zeros({natoms, 3}, torch::kDouble);

    // extract the atomic numbers and positions
    for (int i = 0; i < natoms; i++) {

        // extract the atomic number and coordinates
        double x, y, z; std::getline(system_file_stream, line), std::stringstream(line) >> sm >> x >> y >> z;

        // assign the atomic number and position
        coordinates.index_put_({i, 0}, x), coordinates.index_put_({i, 1}, y), coordinates.index_put_({i, 2}, z), atomic_numbers.at(i) = symbol_to_atomic_number_map(sm) + 1;
    }
}

int System::basis_functions() const {
    return get_shells().nbf();
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

libint2::BasisSet System::get_shells() const {
    return libint2::BasisSet(basis, get_atoms());
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

int System::occupied_spatial_orbitals() const {
    return std::accumulate(atomic_numbers.begin(), atomic_numbers.end(), 0) / 2;
}

int System::occupied_spinorbitals() const {
    return 2 * occupied_spatial_orbitals();
}

int System::virtual_spinorbitals() const {
    return 2 * basis_functions() - occupied_spinorbitals();
}

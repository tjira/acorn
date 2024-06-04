#include "system.h"
#include <fstream>

#define PTABLE "H He Li Be B C N O F Ne"

int System::nocc() const {return std::accumulate(AN.data(), AN.data() + AN.rows(), 0) / 2;}

System::System(const std::string& path) {
    // open the file stream, define the periodic table and initialize container variables
    std::ifstream file(path); std::vector<std::string> ptable; std::string line, sm; int natoms;

    // check if the file stream is good
    if (!file.good()) throw std::runtime_error("SYSTEM FILE `" + path + "` DOES NOT EXIST");

    // fill the periodic table
    for (std::stringstream pss(PTABLE); pss.good(); ptable.push_back(""), pss >> ptable.back()) {}

    // define the symbol to atomic number map
    auto SM2AN = [&ptable](const std::string& sm) {return std::find(ptable.begin(), ptable.end(), sm) - ptable.begin();};

    // extract the number of atoms, skip the first two lines and initialize the atomic number and position matrices
    std::getline(file, line); std::stringstream(line) >> natoms; std::getline(file, line); AN = Matrix(natoms, 1); R = Matrix(natoms, 3);

    // extract the atomic numbers and positions
    for (int i = 0; i < natoms; i++) std::getline(file, line), std::stringstream(line) >> sm >> R(i, 0) >> R(i, 1) >> R(i, 2), AN(i, 0) = SM2AN(sm) + 1;
}

double System::nuclearRepulsion() const {
    // create the variable
    double repulsion = 0;

    // calculate the repulsion
    for (int i = 0; i < AN.rows(); i++) for (int j = 0; j < i; j++) {
        double d = (R.row(i) - R.row(j)).norm(); repulsion += AN(i, 0) * AN(j, 0) / d / 1.889726124626;
    }

    // return the nuclear-nuclear repulsion of the system
    return repulsion;
}

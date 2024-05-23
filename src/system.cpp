#include "ptable.h"
#include "system.h"

#include <fstream>

// variable getters
EigenMatrix<> System::atoms() const {return AN;}
EigenMatrix<> System::coords() const {return R;}

// information getters
int System::nocc() const {return std::accumulate(AN.data(), AN.data() + AN.rows(), 0) / 2;}

System::System(const std::string& path) {
    // open the file stream
    std::ifstream file(path);

    // check if the file stream is good
    if (!file.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // initialize the variables
    int natoms; std::string line, sm;

    // extract the number of atoms and skip the first two lines
    std::getline(file, line); std::stringstream(line) >> natoms; std::getline(file, line);

    // initialize the atomic number and position matrices
    AN = EigenMatrix<>(natoms, 1); R = EigenMatrix<>(natoms, 3);

    // extract the atomic numbers and positions
    for (int i = 0; i < natoms; i++) {
        std::getline(file, line); std::stringstream(line) >> sm >> R(i, 0) >> R(i, 1) >> R(i, 2); AN(i, 0) = SM2AN[sm];
    }
}


double System::nuclearRepulsion() const {
    // create the variable
    double repulsion = 0;

    // calculate the repulsion
    for (int i = 0; i < AN.rows(); i++) {
        for (int j = 0; j < i; j++) {
            double d = (R.row(i) - R.row(j)).norm(); repulsion += AN(i, 0) * AN(j, 0) / d / A2BOHR;
        }
    }

    return repulsion; // return the nuclear-nuclear repulsion of the system
}

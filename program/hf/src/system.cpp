#include "system.h"

#define PTABLE "H He Li Be B C N O F Ne"

int System::nocc() const {return std::accumulate(AN.data(), AN.data() + AN.size(), 0) / 2;}

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
    std::getline(file, line); std::stringstream(line) >> natoms; std::getline(file, line); AN.resize(natoms); R = torch::zeros({natoms, 3}, torch::kDouble);

    // extract the atomic numbers and positions
    for (int i = 0; i < natoms; i++) {
        double x, y, z; std::getline(file, line), std::stringstream(line) >> sm >> x >> y >> z; R.index_put_({i, 0}, x), R.index_put_({i, 1}, y), R.index_put_({i, 2}, z), AN.at(i) = SM2AN(sm) + 1;
    }
}

double System::nuclearRepulsion() const {
    // create the variable
    double repulsion = 0;

    // calculate the repulsion
    for (int i = 0; i < AN.size(); i++) for (int j = 0; j < i; j++) {
        double d = (R.index({i, None}) - R.index({j, None})).norm().item<double>(); repulsion += AN.at(i) * AN.at(j) / d / 1.889726124626;
    }

    // return the nuclear-nuclear repulsion of the system
    return repulsion;
}

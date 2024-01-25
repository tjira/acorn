#include "modelsystem.h"

ModelSystem::ModelSystem(int m, const std::vector<std::vector<std::string>>& potential, const std::vector<std::vector<double>> limits, int ngrid) : limits(limits), ngrid(ngrid), m(m), potential(potential) {}

void ModelSystem::SaveWavefunction(const std::string& fname, const Vector<>& r, const std::vector<Vector<std::complex<double>>>& wfns, const std::vector<double>& energy) {
    // open the file
    std::ofstream file(fname);

    // write the wavefunction
    for (int i = 0; i < wfns.size(); i++) {
        // set the precision and print the independent variable value header
        file << std::fixed << std::setprecision(14) << "#" << std::setw(7) << std::setfill('0') << i
             << " r1                  real                 imag         E=" << energy.at(i) << "\n";

        // write the wavefunction values
        for (long int j = 0; j < r.size(); j++) {
            file << std::setw(20) << std::setfill(' ') << r(j) << " " << std::setw(20) << wfns.at(i)(j).real() << " " << std::setw(20) << wfns.at(i)(j).imag() << "\n";
        }
    }
}

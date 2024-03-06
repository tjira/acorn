#include "modelsystem.h"

ModelSystem::ModelSystem(int m, const std::vector<std::vector<std::string>>& potential, const std::vector<double> limits, int ngrid) : potential(potential), limits(limits), ngrid(ngrid), m(m) {}

void ModelSystem::SaveWavefunction(const std::string& fname, const Vector<>& x, const std::vector<Matrix<std::complex<double>>>& wfns, const std::vector<double>& energy) {
    // open the file
    std::ofstream file(fname);

    // write the wavefunction
    for (size_t i = 0; i < wfns.size(); i++) {
        // set the precision and print the independent variable value header
        file << std::fixed << std::setprecision(14) << "#" << std::setw(7) << std::setfill('0') << i
             << " r1                  real                 imag         E=" << energy.at(i) << "\n";

        // write the wavefunction values
        for (long int j = 0; j < x.size(); j++) {
            file << std::setw(20) << std::setfill(' ') << x(j) << " " << std::setw(20) << wfns.at(i)(j).real() << " " << std::setw(20) << wfns.at(i)(j).imag() << "\n";
        }
    }
}

void ModelSystem::SaveWavefunction(const std::string& fname, const Vector<>& x, const Vector<>& y, const std::vector<Matrix<std::complex<double>>>& wfns, const std::vector<double>& energy) {
    // open the file
    std::ofstream file(fname);

    // write the wavefunction
    for (size_t i = 0; i < wfns.size(); i++) {
        // set the precision and print the independent variable value header
        file << std::fixed << std::setprecision(14) << "#" << std::setw(7) << std::setfill('0') << i
             << " r1                   r2                   real                 imag         E=" << energy.at(i) << "\n";

        // write the wavefunction values
        for (long int j = 0; j < x.size(); j++) {
            for (long int k = 0; k < y.size(); k++) {
                file << std::setw(20) << std::setfill(' ') << x(j) << " " << std::setw(20) << y(k) << std::setw(20) << wfns.at(i)(j, k).real() << " " << std::setw(20) << wfns.at(i)(j, k).imag() << "\n";
            }
        }
    }
}

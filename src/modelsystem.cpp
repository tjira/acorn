#include "modelsystem.h"

ModelSystem::ModelSystem(const std::string& potential, const std::vector<std::vector<double>> limits, int ngrid) : limits(limits), potential(potential), ngrid(ngrid) {}

void ModelSystem::SaveWavefunction(const std::string& fname, const Vector<>& r, const std::vector<std::vector<Vector<std::complex<double>>>>& wfns, const Vector<>& energy) {
    // define the vector size and open the file
    int size = 0; std::ofstream file(fname);

    // calculate the maximum state length 
    for (size_t i = 0; i < wfns.size(); i++) {
        size = std::max(size, (int)wfns.at(i).size());
    }

    // write the wavefunction
    for (int i = 0; i < size; i++) {
        // set the precision and print the independent variable value header
        file << std::fixed << std::setprecision(14) << "#        r1         ";

        // write the header for all states
        for (size_t j = 0; j < wfns.size(); j++) {
            file << "     state" << (j < 10 ? "0" : "") << j << ".real         " << "state" << (j < 10 ? "0" : "") << j << ".imag    ";
        }

        // write the wavefunction energy
        if (i == 0) {
            file << " E=["; for (int j = 0; j < energy.size(); j++) file << energy(j) << (j < energy.size() - 1 ? ", " : "]");
        }

        // write the wavefunction values
        for (long int j = 0; j < r.size(); j++) {
            // print the independent vaiables
            file << (j == 0 ? "\n" : "") << std::setw(20) << r(j);

            // print the dependent variables
            for (size_t k = 0; k < wfns.size(); k++) {
                file << " " << std::setw(20) << wfns.at(k).at(std::min(i, (int)wfns.at(k).size() - 1))(j).real()
                << " " << std::setw(20) << wfns.at(k).at(std::min(i, (int)wfns.at(k).size() - 1))(j).imag();
            }

            // print newline
            file << "\n";
        }
    }
}

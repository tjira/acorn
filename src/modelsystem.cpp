#include "modelsystem.h"

ModelSystem::ModelSystem(int m, const std::vector<std::vector<std::string>>& potential, const std::vector<double> limits, int ngrid) : potential(potential), limits(limits), ngrid(ngrid), m(m) {}

void ModelSystem::SaveWavefunction(const std::filesystem::path& fname, const Vector<>& x, const std::vector<Matrix<std::complex<double>>>& wfns, const std::vector<double>& energy) {
    // open the file, write the header and set the precision
    std::ofstream file(fname); file << "# 1 " << wfns.size() << " " << x.size() << " " << wfns.at(0).size() << std::endl << std::fixed << std::setprecision(12);

    // write the independent coordinates
    for (int i = 0; i < x.size(); i++) {file << (i ? " " : "") << x(i);} file << std::endl;

    // write the wavefunction
    for (size_t i = 0; i < wfns.size(); i++) {
        // the header
        file << "#" << i << " E=" << energy.at(i) << std::endl;

        // the real part
        for (int j = 0; j < wfns.at(i).rows(); j++) {
            for (int k = 0; k < wfns.at(i).cols(); k++) {
                file << (j || k ? " " : "") << wfns.at(i)(j, k).real();
            }
        } file << std::endl;

        // the imaginary part
        for (int j = 0; j < wfns.at(i).rows(); j++) {
            for (int k = 0; k < wfns.at(i).cols(); k++) {
                file << (j || k ? " " : "") << wfns.at(i)(j, k).imag();
            }
        } file << std::endl;
    }
}

void ModelSystem::SaveWavefunction(const std::filesystem::path& fname, const Vector<>& x, const Vector<>& y, const std::vector<Matrix<std::complex<double>>>& wfns, const std::vector<double>& energy) {
    // open the file, write the header and set the precision
    std::ofstream file(fname); file << "# 2 " << wfns.size() << " " << x.size() << " " << y.size() << " " << wfns.at(0).size() << std::endl << std::fixed << std::setprecision(12);

    // write the independent coordinates
    for (int i = 0; i < x.size(); i++) {file << (i ? " " : "") << x(i);} file << std::endl;
    for (int i = 0; i < y.size(); i++) {file << (i ? " " : "") << y(i);} file << std::endl;

    // write the wavefunction
    for (size_t i = 0; i < wfns.size(); i++) {
        // the header
        file << "#" << i << " E=" << energy.at(i) << std::endl;

        // the real part
        for (int j = 0; j < wfns.at(i).rows(); j++) {
            for (int k = 0; k < wfns.at(i).cols(); k++) {
                file << (j || k ? " " : "") << wfns.at(i)(j, k).real();
            }
        } file << std::endl;

        // the imaginary part
        for (int j = 0; j < wfns.at(i).rows(); j++) {
            for (int k = 0; k < wfns.at(i).cols(); k++) {
                file << (j || k ? " " : "") << wfns.at(i)(j, k).imag();
            }
        } file << std::endl;
    }
}

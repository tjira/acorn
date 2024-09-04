#pragma once

#include       "expression.h"
#include "fouriertransform.h"
#include            "input.h"

class Wavefunction {
public:
    Wavefunction() = default; Wavefunction(const Input::Wavefunction& input);

    Wavefunction operator-(const Wavefunction& other_wavefunction) const;
    Wavefunction operator*(const std::complex<double>& scalar) const;
    Wavefunction adiabatized(const std::vector<Eigen::MatrixXd>& transformation_matrices) const;
    double energy(const Eigen::MatrixXd& diabatic_potential, const Eigen::MatrixXd& fourier_grid) const;
    Eigen::VectorXd momentum(const Eigen::MatrixXd& fourier_grid) const;
    Eigen::VectorXd position(const Eigen::MatrixXd& grid) const;
    Eigen::MatrixXcd get_data() const;
    Eigen::MatrixXd get_density() const;
    Eigen::MatrixXd get_grid() const;
    double get_grid_spacing() const;
    Eigen::MatrixXd get_fourier_grid() const;
    std::vector<int> get_shape() const;
    std::vector<std::string> get_variables() const;
    Wavefunction get_normalized() const;
    std::complex<double> overlap(const Wavefunction& other_wavefunction) const;
    Wavefunction propagated(const std::vector<Eigen::MatrixXcd>& real_propagators, const std::vector<Eigen::MatrixXcd>& fourier_propagators) const;
    std::tuple<std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>> propagators(const Eigen::MatrixXd& diabatic_potential, const Eigen::MatrixXd& fourier_grid, const std::complex<double>& unit, double step) const;

private:
    Input::Wavefunction input; Eigen::MatrixXcd data;
};

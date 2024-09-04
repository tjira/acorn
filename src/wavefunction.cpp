#include "wavefunction.h"

std::tuple<Wavefunction, Eigen::MatrixXd> Wavefunction::Initialize(const Input::Wavefunction& input) {
    // initialize the wavefunction
    Wavefunction wavefunction; wavefunction.input = input, wavefunction.data = Eigen::MatrixXd((int)std::pow(input.grid_points, input.dimension), input.guess.size() / 2);

    // calculate the grid and get the variables
    Eigen::MatrixXd grid = wavefunction.get_grid(); std::vector<std::string> variables = wavefunction.get_variables();

    // evaluate the guess at the grid points
    for (int i = 0; i < wavefunction.data.cols(); i++) {
        wavefunction.data.col(i) = Expression(input.guess.at(2 * i), variables).evaluate(grid) + std::complex<double>(0, 1) * Expression(input.guess.at(2 * i + 1), variables).evaluate(grid);
    }

    // add the momentum to the wavefunction
    for (int i = 0; i < wavefunction.data.cols(); i++) {
        wavefunction.data.col(i) = wavefunction.data.col(i).array() * (std::complex<double>(0, 1) * input.momentum * grid.rowwise().sum().array()).exp();
    }

    // normalize the wavefunction
    wavefunction.data /= std::sqrt(std::abs(wavefunction.overlap(wavefunction)));

    // return the wavefunction and the grid
    return {wavefunction, grid};
}

Wavefunction Wavefunction::operator*(const std::complex<double>& scalar) const {
    Wavefunction multiplied_wavefunction(*this); multiplied_wavefunction.data *= scalar; return multiplied_wavefunction;
}

Wavefunction Wavefunction::operator-(const Wavefunction& other_wavefunction) const {
    Wavefunction subtracted_wavefunction(*this); subtracted_wavefunction.data -= other_wavefunction.data; return subtracted_wavefunction;
}

Wavefunction Wavefunction::adiabatized(const std::vector<Eigen::MatrixXd>& transformation_matrices) const {
    // define adiabatic wavefunction
    Wavefunction adiabatized_wavefunction(*this);

    // loop over all points and transform the wavefunction
    for (int j = 0; j < data.rows(); j++) adiabatized_wavefunction.data.row(j) = transformation_matrices.at(j).adjoint() * data.row(j).transpose();

    // return the wfn
    return adiabatized_wavefunction;
}

double Wavefunction::energy(const Eigen::MatrixXd& diabatic_potential, const Eigen::MatrixXd& fourier_grid) const {
    // define the kinetic and potential energy
    double Ek = 0, Ep = 0;

    // loop over all wavefunctions
    for (int i = 0, n = data.cols(); i < data.cols(); i++) {

        // calculate the kinetic energy contribution
        Ek += (data.col(i).adjoint() * FourierTransform::IFFT(fourier_grid.array().pow(2).rowwise().sum() * FourierTransform::FFT(data.col(i), get_shape()).array(), get_shape()))(0).real();

        // calculate the potential energy contributions
        for (int j = 0; j < n; j++) Ep += (data.col(i).adjoint() * (diabatic_potential.col(i * n + j).array() * data.col(j).array()).matrix())(0).real();
    }

    // return the total energy
    return (0.5 * Ek / input.mass + Ep) * get_grid_spacing();
}

Eigen::VectorXd Wavefunction::momentum(const Eigen::MatrixXd& fourier_grid) const {
    // define the momentum container
    Eigen::VectorXd momentum = Eigen::VectorXd::Zero(input.dimension);

    // calculate the momentum
    for (int i = 0; i < data.cols(); i++) {
        momentum.array() -= (data.col(i).replicate(1, input.dimension).conjugate().array() * FourierTransform::IFFT(fourier_grid.rowwise().sum().array() * FourierTransform::FFT(data.col(i), get_shape()).array(), get_shape()).replicate(1, input.dimension).array()).sum().real();
    }

    // return the momentum
    return momentum * get_grid_spacing();
}

Eigen::VectorXd Wavefunction::position(const Eigen::MatrixXd& grid) const {
    // define the position container
    Eigen::VectorXd position = Eigen::VectorXd::Zero(input.dimension);

    // calculate the position
    for (int i = 0; i < data.cols(); i++) {
        position.array() += (data.col(i).replicate(1, input.dimension).conjugate().array() * grid.array() * data.col(i).replicate(1, input.dimension).array()).colwise().sum().real();
    }

    // return the position
    return position * get_grid_spacing();
}

Eigen::MatrixXcd Wavefunction::get_data() const {
    return data;
}

Eigen::MatrixXd Wavefunction::get_density() const {
    return (data.transpose() * data.conjugate()).array().abs() * get_grid_spacing();
}

Eigen::MatrixXd Wavefunction::get_grid() const {
    // create the container for the grid
    Eigen::MatrixXd grid((int)std::pow(input.grid_points, input.dimension), input.dimension);

    // fill the grid
    for (int i = 0; i < grid.rows(); i++) for (int j = 0; j < grid.cols(); j++) {
        grid(i, j) = input.grid_limits.at(0) + ((i / (int)std::pow(input.grid_points, input.dimension - j - 1)) % input.grid_points) * get_grid_spacing();
    }

    // return the grid
    return grid;
}

double Wavefunction::get_grid_spacing() const {
    return (input.grid_limits.at(1) - input.grid_limits.at(0)) / (input.grid_points - 1);
}

Eigen::MatrixXd Wavefunction::get_fourier_grid() const {
    // create the container for the grid
    Eigen::MatrixXd fourier_grid((int)std::pow(input.grid_points, input.dimension), input.dimension);

    // fill the constant contributions to the last column of the fourier grid 
    fourier_grid.col(input.dimension - 1).fill(2 * M_PI / input.grid_points / get_grid_spacing());

    // fill the non-constant contributions to the last column of the fourier grid
    for (int i = 0; i < fourier_grid.rows() / input.grid_points; i++) for (int j = 0; j < input.grid_points; j++) {
        fourier_grid.col(input.dimension - 1)(i * input.grid_points + j) *= j - (j < input.grid_points / 2 ? 0 : input.grid_points);
    }

    // fill the rest of the columns of the fourier grid
    for (int i = 0; i < fourier_grid.rows(); i++) for (int j = 0; j < fourier_grid.cols(); j++) {
        fourier_grid(i, j) = fourier_grid(i / (int)std::pow(input.grid_points, input.dimension - j - 1), input.dimension - 1);
    }

    // return the grid
    return fourier_grid;
}

Wavefunction Wavefunction::get_normalized() const {
    Wavefunction wfn(*this); wfn.data /= std::sqrt(std::abs(overlap(wfn))); return wfn;
}

std::vector<int> Wavefunction::get_shape() const {
    return std::vector<int>(input.dimension, input.grid_points);
}

std::vector<std::string> Wavefunction::get_variables() const {
    // create the container for the variables
    std::vector<std::string> variables(input.dimension);

    // fill the variables
    for (int i = 0; i < input.dimension; i++) variables.at(i) = input.dimension > 3 ? "r" + std::to_string(i + 1) : std::vector<std::string>{"x", "y", "z"}.at(i);

    // return the variables
    return variables;
}

std::complex<double> Wavefunction::overlap(const Wavefunction& other_wavefunction) const {
    // define the overlap container
    std::complex<double> overlap = 0;
    
    // calculate all the contributions to the overlap
    for (int i = 0; i < data.cols(); i++) for (int j = 0; j < other_wavefunction.data.cols(); j++) overlap += (data.col(i).adjoint() * other_wavefunction.data.col(j))(0) * get_grid_spacing();

    // return the overlap
    return overlap;
}

Wavefunction Wavefunction::propagated(const std::vector<Eigen::MatrixXcd>& R, const std::vector<Eigen::MatrixXcd>& K) const {
    // output wavefunction
    Wavefunction propagated_wavefunction(*this);

    // perform the half step in the real space
    for (int i = 0; i < data.rows(); i++) propagated_wavefunction.data.row(i) = R.at(i) * propagated_wavefunction.data.row(i).transpose();

    // perform the fourier transform
    for (int j = 0; j < data.cols(); j++) propagated_wavefunction.data.col(j) = FourierTransform::FFT(propagated_wavefunction.data.col(j), get_shape());

    // perform the full step in the momentum space
    for (int i = 0; i < data.rows(); i++) propagated_wavefunction.data.row(i) = K.at(i) * propagated_wavefunction.data.row(i).transpose();

    // perform the inverse fourier transform
    for (int j = 0; j < data.cols(); j++) propagated_wavefunction.data.col(j) = FourierTransform::IFFT(propagated_wavefunction.data.col(j), get_shape());

    // perform the half step in the real space
    for (int i = 0; i < data.rows(); i++) propagated_wavefunction.data.row(i) = R.at(i) * propagated_wavefunction.data.row(i).transpose();

    // return the propagated wavefuction
    return propagated_wavefunction;
}

std::tuple<std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>> Wavefunction::propagators(const Eigen::MatrixXd& diabatic_potential, const Eigen::MatrixXd& fourier_grid, const std::complex<double>& unit, double step) const {
    // define the propagator array
    std::vector<Eigen::MatrixXcd> real_propagators(diabatic_potential.rows()), fourier_propagators(diabatic_potential.rows());

    // calculate the propagators
    for (int i = 0, n = data.cols(); i < diabatic_potential.rows(); i++) {

        // initialize the eigenvalue solver and diagonalize the potential
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(diabatic_potential.row(i).reshaped(n, n)); Eigen::MatrixXcd C = solver.eigenvectors(), E = solver.eigenvalues();

        // momentum space propagator as a diagonal matrix with kinetic energy operators
        Eigen::MatrixXcd fourier_propagator = (-0.5 * unit * fourier_grid.row(i).array().pow(2).sum() * step * Eigen::VectorXd::Ones(n).array() / input.mass).exp();

        // real space propagator as a general matrix exponential
        Eigen::MatrixXcd real_propagator = C * (-0.5 * unit * E.array() * step).exp().matrix().asDiagonal() * C.adjoint();

        // store the propagator matrices
        real_propagators.at(i) = real_propagator, fourier_propagators.at(i) = fourier_propagator.asDiagonal();
    }

    // return the propagator
    return {real_propagators, fourier_propagators};
}

#include "fourier.h"
#include "wavefunction.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Nonadiabatic Quantum Dynamics Program", "1.0", argparse::default_arguments::none);

    // add the command line arguments
    program.add_argument("-g", "--guess").help("-- Initial wavefunction.").nargs(2).default_value(std::vector<std::string>{"exp(-(x+10)^2)", "0"});
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Maximum number of iterations.").default_value(350).scan<'i', int>();
    program.add_argument("-m", "--mass").help("-- Mass of the model system.").default_value(2000.0).scan<'g', double>();
    program.add_argument("-p", "--momentum").help("-- Momentum of the model system.").default_value(10.95).scan<'g', double>();
    program.add_argument("-s", "--step").help("-- Time step of the simulation.").default_value(10.0).scan<'g', double>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    int iters = program.get<int>("-i"); double mass = program.get<double>("-m"), step = program.get<double>("-s");

    // load the potential function and define wfn vector
    EigenMatrix<> U = Eigen::LoadMatrix("U.mat"); std::vector<Wavefunction> wfn;

    // initialize the wavefunctions
    wfn.emplace_back(program.get<std::vector<std::string>>("-g").at(0), U.leftCols(1), mass, program.get<double>("-p")); wfn.at(0) = wfn.at(0).norm() ? wfn.at(0) / wfn.at(0).norm() : wfn.at(0);
    wfn.emplace_back(program.get<std::vector<std::string>>("-g").at(1), U.leftCols(1), mass, program.get<double>("-p")); wfn.at(1) = wfn.at(1).norm() ? wfn.at(1) / wfn.at(1).norm() : wfn.at(1);

    // remove first column from matrix
    U.block(0, 0, U.rows(), U.cols() - 1) = U.rightCols(U.cols() - 1); U.conservativeResize(U.rows(), U.cols() - 1);

    // define the real time kinetic and potential operators
    std::vector<std::vector<EigenMatrix<std::complex<double>>>> R(wfn.size(), std::vector<EigenMatrix<std::complex<double>>>(wfn.size()));
    std::vector<std::vector<EigenMatrix<std::complex<double>>>> K(wfn.size(), std::vector<EigenMatrix<std::complex<double>>>(wfn.size()));

    K.at(0).at(0) = (-0.5 * std::complex<double>(0, 1) * wfn.at(0).getk().array().pow(2) * step / mass).exp();
    K.at(1).at(1) = (-0.5 * std::complex<double>(0, 1) * wfn.at(1).getk().array().pow(2) * step / mass).exp();
    EigenVector<std::complex<double>> D = 4 * U.col(2).array().abs().pow(2) + (U.col(0) - U.col(3)).array().pow(2);
    EigenVector<std::complex<double>> a = (-0.25 * std::complex<double>(0, 1) * (U.col(0) + U.col(3)) * step).array().exp();
    EigenVector<std::complex<double>> b = (0.25 * D.array().sqrt() * step).cos();
    EigenVector<std::complex<double>> c = std::complex<double>(0, 1) * (0.25 * D.array().sqrt() * step).sin() / D.array().sqrt();
    R.at(0).at(0) = a.array() * (b.array() + c.array() * (U.col(3) - U.col(0)).array());
    R.at(0).at(1) = -2 * a.array() * c.array() * U.col(1).array();
    R.at(1).at(0) = -2 * a.array() * c.array() * U.col(1).array();
    R.at(1).at(1) = a.array() * (b.array() + c.array() * (U.col(0) - U.col(3)).array());

    double E = 0;

    for (int i = 0; i < iters; i++) {
        // save the previous values and temporary wavefunction
        std::vector<Wavefunction> wfnp = wfn, wfnt = wfn; double Ep = E;

        // apply the propagator
        wfnt.at(0).get() = R.at(0).at(0).array() * wfn.at(0).get().array() + R.at(0).at(1).array() * wfn.at(1).get().array();
        wfnt.at(1).get() = R.at(1).at(0).array() * wfn.at(0).get().array() + R.at(1).at(1).array() * wfn.at(1).get().array(); wfn = wfnt;
        wfnt.at(0).get() = FourierTransform::Inverse(K.at(0).at(0).array() * FourierTransform::Forward(wfn.at(0).get()).array());
        wfnt.at(1).get() = FourierTransform::Inverse(K.at(1).at(1).array() * FourierTransform::Forward(wfn.at(1).get()).array()); wfn = wfnt;
        wfnt.at(0).get() = R.at(0).at(0).array() * wfn.at(0).get().array() + R.at(0).at(1).array() * wfn.at(1).get().array();
        wfnt.at(1).get() = R.at(1).at(0).array() * wfn.at(0).get().array() + R.at(1).at(1).array() * wfn.at(1).get().array(); wfn = wfnt;

        // copy the wavefunction and transform it to the adiabatic basis if the matrices UT were calculated
        // psitemp = psi; for (int j = 0; j < psitemp.rows(); j++) psitemp.row(j) = UT.at(j).transpose() * psitemp.row(j).transpose();

        // calculate the density matrix for the current time step in the chosen basis
        // rho = psitemp.transpose() * psitemp.conjugate() * dr;

        // calculate the total energy
        EigenMatrix<std::complex<double>> Ek00 = Eigen::ComplexConjugate(wfn.at(0).get()).array() * FourierTransform::Inverse(wfn.at(0).getk().array().pow(2) * FourierTransform::Forward(wfn.at(0).get()).array()).array();
        EigenMatrix<std::complex<double>> Ek11 = Eigen::ComplexConjugate(wfn.at(1).get()).array() * FourierTransform::Inverse(wfn.at(1).getk().array().pow(2) * FourierTransform::Forward(wfn.at(1).get()).array()).array();
        EigenMatrix<std::complex<double>> Ep00 = Eigen::ComplexConjugate(wfn.at(0).get()).array() * U.col(0).array() * wfn.at(0).get().array();
        EigenMatrix<std::complex<double>> Ep01 = Eigen::ComplexConjugate(wfn.at(0).get()).array() * U.col(1).array() * wfn.at(1).get().array();
        EigenMatrix<std::complex<double>> Ep10 = Eigen::ComplexConjugate(wfn.at(1).get()).array() * U.col(2).array() * wfn.at(0).get().array();
        EigenMatrix<std::complex<double>> Ep11 = Eigen::ComplexConjugate(wfn.at(1).get()).array() * U.col(3).array() * wfn.at(1).get().array();
        double E = (0.5 * (Ek00 + Ek11) / mass + Ep00 + Ep01 + Ep10 + Ep11).sum().real() * (wfn.at(0).getr()(1) - wfn.at(0).getr()(0));
        // E = Propagator<2>::Energy(system, V.real(), res.msv.r, res.msv.k.array().pow(2), psi);

        // calculate the errors and assign the previous values
        // double Eerr = std::abs(E - Ep), Derr = (psi.array().abs2() - psip.array().abs2()).abs2().sum();

        // append the temporary wfn in adiabatic basis and energy
        // energies.push_back(E); if (optn.savewfn) psis.push_back(psitemp);

        // print the iteration
        std::printf("%6d %20.14f\n", i + 1, E);
    }
    // EigenVector<> a(4); a << 0, 1, 2, 3;

    // reshape to matrix 2x2
    // EigenMatrix<> b = Eigen::Map<EigenMatrix<>>(a.data(), 2, 2);
}

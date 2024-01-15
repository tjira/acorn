#include "method.h"

Vector<> Method::frequency(const System& system, const Matrix<>& H) {
    // create the mass matrix
    Matrix<> MM(3 * system.getAtoms().size(), 3 * system.getAtoms().size());
    for (int i = 0; i < MM.rows(); i++) {
        MM(i, i) = std::sqrt(1 / masses.at(system.getAtoms().at(i / 3).atomic_number));
    }

    // calculate the frequencies in atomic units, convert them, exract them and return
    Eigen::EigenSolver<Matrix<>> solver(MM * H * MM); auto eval = solver.eigenvalues();
    Vector<> freqs = (eval.cwiseSqrt() - eval.cwiseSqrt().imag()).real() * CFREQ;
    std::sort(freqs.begin(), freqs.end(), std::greater<>()); return freqs;
}

Method::Result Method::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // define the gradient matrix
    res.G = Matrix<>(system.getAtoms().size(), 3);

    std::cout << opt.gradient.step << std::endl;

    // print the header
    if (print) std::printf("\n  ELEM      dE [Eh/Bohr]        TIME\n");

    #pragma omp parallel for num_threads(nthread) collapse(2)
    for (int i = 0; i < res.G.rows(); i++) {
        for (int j = 0; j < res.G.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices and temporary systems
            Matrix<> dirMinus(system.getAtoms().size(), 3); System sysMinus = system;
            Matrix<> dirPlus(system.getAtoms().size(), 3); System sysPlus = system;

            // fill the direction matrices
            dirMinus(i, j) -= opt.gradient.step * A2BOHR; dirPlus(i, j) += opt.gradient.step * A2BOHR;

            // move the systems
            sysMinus.move(dirMinus), sysPlus.move(dirPlus);

            // calculate and assign the derivative
            res.G(i, j) = BOHR2A * (energy(sysPlus, res) - energy(sysMinus, res)) / opt.gradient.step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // return the gradient
    return res;
}

Method::Result Method::hessian(const System& system, const Integrals&, Result res, bool print) const {
    // define the gradient matrix
    res.H = Matrix<>(3 * system.getAtoms().size(), 3 * system.getAtoms().size());

    // print the header
    if (print) std::printf("\n  ELEM      dE [Eh/Bohr]        TIME\n");

    #pragma omp parallel for num_threads(nthread) collapse(2)
    for (int i = 0; i < res.H.rows(); i++) {
        for (int j = 0; j < res.H.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices and temporary systems
            System sysMinusMinus = system, sysMinusPlus = system, sysPlusMinus = system, sysPlusPlus = system;
            Matrix<> dir1(system.getAtoms().size(), 3), dir2(system.getAtoms().size(), 3);

            // fill the direction matrices
            dir1(i / 3, i % 3) = opt.hessian.step * A2BOHR; dir2(j / 3, j % 3) = opt.hessian.step * A2BOHR;

            // move the systems
            sysMinusMinus.move(-dir1 - dir2), sysMinusPlus.move(-dir1 + dir2);
            sysPlusMinus.move(dir1 - dir2), sysPlusPlus.move(dir1 + dir2);

            // calculate the energies
            double energyMinusMinus = energy(sysMinusMinus, res), energyMinusPlus = energy(sysMinusPlus, res);
            double energyPlusMinus = energy(sysPlusMinus, res), energyPlusPlus = energy(sysPlusPlus, res);

            // calculate and assign the derivative
            res.H(i, j) = BOHR2A * BOHR2A * (energyPlusPlus - energyMinusPlus - energyPlusMinus + energyMinusMinus) / opt.hessian.step / opt.hessian.step / 4;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // return the gradient
    return res;
}

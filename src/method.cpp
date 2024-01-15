#include "method.h"

Method::Result Method::gradient(const System& system, const Integrals& ints, Result res, bool print) const {
    // define the gradient matrix
    res.G = Matrix<>(system.getAtoms().size(), 3);

    double step = 1e-6;

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
            dirMinus(i, j) -= step * A2BOHR; dirPlus(i, j) += step * A2BOHR;

            // move the systems
            sysMinus.move(dirMinus), sysPlus.move(dirPlus);

            // calculate and assign the derivative
            res.G(i, j) = BOHR2A * (energy(sysPlus) - energy(sysMinus)) / step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // return the gradient
    return res;
}

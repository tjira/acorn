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

            // move the systems and create the integral containers
            sysMinus.move(dirMinus), sysPlus.move(dirPlus); Integrals intsMinus, intsPlus;

            // calculate all the atomic integrals
            libint2::initialize();
            intsMinus.S = Integral::Overlap(sysMinus), intsMinus.T = Integral::Kinetic(sysMinus);
            intsMinus.V = Integral::Nuclear(sysMinus), intsMinus.J = Integral::Coulomb(sysMinus);
            intsPlus.S = Integral::Overlap(sysPlus), intsPlus.T = Integral::Kinetic(sysPlus);
            intsPlus.V = Integral::Nuclear(sysPlus), intsPlus.J = Integral::Coulomb(sysPlus);
            libint2::finalize();

            // calculate and assign the derivative
            res.G(i, j) = BOHR2A * (run(sysPlus, intsPlus, {}, false).Etot - run(sysMinus, intsMinus, {}, false).Etot) / step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // return the gradient
    return res;
}

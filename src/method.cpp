#include "restrictedmollerplesset.h"

// include interfaces
#include "orca.h"

template <class M>
void Method<M>::dynamics(System system, const Integrals& ints, Result res, bool print) const {
    // print the header
    if (print) std::printf("\n ITER  TIME [fs]        E [Eh]              KIN [Eh]              TK [K]         |GRAD|      TIME\n");

    // create position, velocity and mass matrices
    Matrix<> v(system.getAtoms().size(), 3), a(system.getAtoms().size(), 3), m(system.getAtoms().size(), 3);

    // fill the mass matrix
    for(size_t i = 0; i < system.getAtoms().size(); i++) {
        m.row(i) = [](double m) {return Vector<>::Constant(3, m);}(an2m.at(system.getAtoms().at(i).atomic_number));
    }

    // get the degrees of freedom and write the initial geometry
    double Nf = system.getAtoms().size() * 3 - 6; system.save(get()->opt.dynamics.output) ;

    // start the timer
    auto start = Timer::Now();

    // calculate the initial energy with gradient
    if constexpr (std::is_same_v<M, Orca>) {
        res = gradient(system, ints, res, false);
    } else {
        res = run(system, res, false); res = gradient(system, ints, res, false);
    }

    // calculate the initial kinetic energy and temperature
    double Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(); double T = 2 * Ekin / Nf;

    // print the zeroth iteration
    if (print) std::printf("%6d %9.4f %20.14f %20.14f %20.14f %.2e %s\n", 0, 0.0, res.Etot, Ekin, T / BOLTZMANN, res.G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());

    for (int i = 0; i < get()->opt.dynamics.iters; i++) {
        // start the timer and store the previous v and a
        start = Timer::Now(); Matrix<> vp = v, ap = a;

        // calculate the velocity and accceleration
        a = -res.G.array() / m.array(); v = vp + (ap + a) * get()->opt.dynamics.step / 2;

        // calculate the kinetic energy and temperature
        Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(); T = 2 * Ekin / Nf;

        // move the system
        system.move(get()->opt.dynamics.step * (v + 0.5 * a * get()->opt.dynamics.step));

        // write the current geometry
        system.save(get()->opt.dynamics.output, std::ios::app);

        // calculate the next energy with gradient
        if constexpr (std::is_same_v<M, Orca>) {
            res = gradient(system, ints, res, false);
        } else {
            res = run(system, res, false); res = gradient(system, ints, res, false);
        }

        // print the iteration info
        if (print) std::printf("%6d %9.4f %20.14f %20.14f %20.14f %.2e %s\n", i + 1, AU2FS * get()->opt.dynamics.step * (i + 1), res.Etot, Ekin, T / BOLTZMANN, res.G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

}

template <class M>
Vector<> Method<M>::frequency(const System& system, const Matrix<>& H) {
    // create the mass matrix
    Matrix<> MM(3 * system.getAtoms().size(), 3 * system.getAtoms().size());
    for (int i = 0; i < MM.rows(); i++) {
        MM(i, i) = std::sqrt(1 / an2m.at(system.getAtoms().at(i / 3).atomic_number));
    }

    // calculate the frequencies in atomic units, convert them, exract them and return
    Eigen::EigenSolver<Matrix<>> solver(MM * H * MM); auto eval = solver.eigenvalues();
    Vector<> freqs = (eval.cwiseSqrt() - eval.cwiseSqrt().imag()).real() * CFREQ;
    std::sort(freqs.begin(), freqs.end(), std::greater<>()); return freqs;
}

template <class M>
Result Method<M>::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // return if run from interface
    if constexpr (std::is_same_v<M, Orca>) return res;

    // define the gradient matrix
    res.G = Matrix<>(system.getAtoms().size(), 3);

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
            dirMinus(i, j) -= get()->opt.gradient.step * A2BOHR; dirPlus(i, j) += get()->opt.gradient.step * A2BOHR;

            // move the systems
            sysMinus.move(dirMinus), sysPlus.move(dirPlus);

            // calculate and assign the derivative
            res.G(i, j) = BOHR2A * (run(sysPlus, res, false).Etot - run(sysMinus, res, false).Etot) / get()->opt.gradient.step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // assign the gradient to the correct method variable
    if constexpr (std::is_same_v<M, RestrictedMollerPlesset>) res.rmp.G = res.G;
    if constexpr (std::is_same_v<M, RestrictedHartreeFock>) res.rhf.G = res.G;

    // return the gradient
    return res;
}

template <class M>
Result Method<M>::hessian(const System& system, const Integrals&, Result res, bool print) const {
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
            dir1(i / 3, i % 3) = get()->opt.hessian.step * A2BOHR; dir2(j / 3, j % 3) = get()->opt.hessian.step * A2BOHR;

            // move the systems
            sysMinusMinus.move(-dir1 - dir2), sysMinusPlus.move(-dir1 + dir2);
            sysPlusMinus.move(dir1 - dir2), sysPlusPlus.move(dir1 + dir2);

            // calculate the energies
            double energyMinusMinus = run(sysMinusMinus, res, false).Etot, energyMinusPlus = run(sysMinusPlus, res, false).Etot;
            double energyPlusMinus = run(sysPlusMinus, res, false).Etot, energyPlusPlus = run(sysPlusPlus, res, false).Etot;

            // calculate and assign the derivative
            res.H(i, j) = BOHR2A * BOHR2A * (energyPlusPlus - energyMinusPlus - energyPlusMinus + energyMinusMinus) / get()->opt.hessian.step / get()->opt.hessian.step / 4;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }


    // assign the hessian to the correct method variable
    if constexpr (std::is_same_v<M, RestrictedMollerPlesset>) res.rmp.H = res.H;
    if constexpr (std::is_same_v<M, RestrictedHartreeFock>) res.rhf.H = res.H;

    // return the gradient
    return res;
}

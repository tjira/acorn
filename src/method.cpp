#include "restrictedconfigurationinteraction.h"

// include interfaces
#include "orca.h"

template <class M>
void Method<M>::dynamics(System system, Integrals ints, Result res, bool print) const {
    // print the header
    if (print) std::printf(" ITER  TIME [fs]    POT [Eh]       KIN [Eh]        E [Eh]         TK [K]      |GRAD|      TIME\n");

    // create position, velocity and mass matrices
    Matrix<> v(system.getAtoms().size(), 3), a(system.getAtoms().size(), 3), m(system.getAtoms().size(), 3);

    // fill the mass matrix
    for(size_t i = 0; i < system.getAtoms().size(); i++) {
        m.row(i) = [](double m) {return Vector<>::Constant(3, m);}(an2m.at(system.getAtoms().at(i).atomic_number) * AMU);
    }

    // get the degrees of freedom and write the initial geometry
    double Nf = system.getAtoms().size() * 3 - 6; system.save(get()->opt.dynamics.folder / std::filesystem::path("trajectory.xyz")) ;

    // get the dynamics step and start the timer
    double step = get()->opt.dynamics.step; auto start = Timer::Now();

    // calculate the initial energy with gradient
    if constexpr (std::is_same_v<M, Orca>) {
        res = gradient(system, ints, res, false);
    } else if constexpr (std::is_same_v<M, RestrictedHartreeFock>) {
        std::tie(res, ints) = run(system, res, false);
        ints.dS = Integral::dOverlap(system), ints.dT = Integral::dKinetic(system);
        ints.dV = Integral::dNuclear(system), ints.dJ = Integral::dCoulomb(system);
        res = gradient(system, ints, res, false);
    } else {
        std::tie(res, ints) = run(system, res, false); res = gradient(system, ints, res, false);
    }

    // calculate the initial kinetic energy and temperature
    double Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(); double T = 2 * Ekin / Nf;

    // print the zeroth iteration
    if (print) std::printf("%6d %9.4f %14.8f %14.8f %14.8f %14.8f %.2e %s\n", 0, 0.0, res.Etot, Ekin, res.Etot + Ekin, T * AU2K, res.G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());

    for (int i = 0; i < get()->opt.dynamics.iters; i++) {
        // start the timer and store the previous v and a
        start = Timer::Now(); Matrix<> vp = v, ap = a;

        // calculate the velocity and accceleration
        a = -res.G.array() / m.array(); v = vp + 0.5 * (ap + a) * step;

        // calculate the kinetic energy and temperature
        Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(); T = 2 * Ekin / Nf;

        // move the system
        system.move(step * (v + 0.5 * a * step));

        // write the current geometry
        system.save(get()->opt.dynamics.folder / std::filesystem::path("trajectory.xyz"), std::ios::app);

        // calculate the next energy with gradient
        if constexpr (std::is_same_v<M, Orca>) {
            res = gradient(system, ints, res, false);
        } else if constexpr (std::is_same_v<M, RestrictedHartreeFock>) {
            std::tie(res, ints) = run(system, res, false);
            ints.dS = Integral::dOverlap(system), ints.dT = Integral::dKinetic(system);
            ints.dV = Integral::dNuclear(system), ints.dJ = Integral::dCoulomb(system);
            res = gradient(system, ints, res, false);
        } else {
            std::tie(res, ints) = run(system, res, false); res = gradient(system, ints, res, false);
        }

        // print the iteration info
        if (print) std::printf("%6d %9.4f %14.8f %14.8f %14.8f %14.8f %.2e %s\n", i + 1, AU2FS * step * (i + 1), res.Etot, Ekin, res.Etot + Ekin, T * AU2K, res.G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
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
    // define the gradient matrix
    res.G = Matrix<>(system.getAtoms().size(), 3);

    // define the step
    double step = 0; if constexpr (!std::is_same_v<M, Orca>) step = get()->opt.gradient.step;

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) collapse(2)
    #endif
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
            res.G(i, j) = BOHR2A * (std::get<0>(run(sysPlus, res, false)).Etot - std::get<0>(run(sysMinus, res, false)).Etot) / step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // assign the gradient to the correct method variable
    if constexpr (std::is_same_v<M, RestrictedConfigurationInteraction>) res.rci.G = res.G;
    if constexpr (std::is_same_v<M, RestrictedMollerPlesset>) res.rmp.G = res.G;
    if constexpr (std::is_same_v<M, RestrictedHartreeFock>) res.rhf.G = res.G;

    // return the gradient
    if (print) {std::cout << std::endl;} return res;
}

template <class M>
Result Method<M>::hessian(const System& system, const Integrals&, Result res, bool print) const {
    // define the gradient matrix
    res.H = Matrix<>(3 * system.getAtoms().size(), 3 * system.getAtoms().size());

    // define the step
    double step = 0; if constexpr (!std::is_same_v<M, Orca>) step = get()->opt.hessian.step;

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < res.H.rows(); i++) {
        for (int j = i; j < res.H.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices and temporary systems
            System sysMinusMinus = system, sysMinusPlus = system, sysPlusMinus = system, sysPlusPlus = system;
            Matrix<> dir1(system.getAtoms().size(), 3), dir2(system.getAtoms().size(), 3);

            // fill the direction matrices
            dir1(i / 3, i % 3) = step * A2BOHR; dir2(j / 3, j % 3) = step * A2BOHR;

            // move the systems
            sysMinusMinus.move(-dir1 - dir2), sysMinusPlus.move(-dir1 + dir2);
            sysPlusMinus.move(dir1 - dir2), sysPlusPlus.move(dir1 + dir2);

            // calculate the energies
            double energyMinusMinus = std::get<0>(run(sysMinusMinus, res, false)).Etot, energyMinusPlus = std::get<0>(run(sysMinusPlus, res, false)).Etot;
            double energyPlusMinus = std::get<0>(run(sysPlusMinus, res, false)).Etot, energyPlusPlus = std::get<0>(run(sysPlusPlus, res, false)).Etot;

            // calculate and assign the derivative
            res.H(i, j) = BOHR2A * BOHR2A * (energyPlusPlus - energyMinusPlus - energyPlusMinus + energyMinusMinus) / step / step / 4; res.H(j, i) = res.H(i, j);

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, res.H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // assign the hessian to the correct method variable
    if constexpr (std::is_same_v<M, RestrictedConfigurationInteraction>) res.rci.H = res.H;
    if constexpr (std::is_same_v<M, RestrictedMollerPlesset>) res.rmp.H = res.H;
    if constexpr (std::is_same_v<M, RestrictedHartreeFock>) res.rhf.H = res.H;

    // print the new line return the gradient
    if (print) {std::cout << std::endl;} return res;
}

// method definitions
template class Method<RestrictedConfigurationInteraction>;
template class Method<RestrictedMollerPlesset>;
template class Method<RestrictedHartreeFock>;

// interface definitions
template class Method<Orca>;

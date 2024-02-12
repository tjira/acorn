#include "restrictedconfigurationinteraction.h"

std::tuple<Result, Integrals> RestrictedConfigurationInteraction::run(const System& system, Result res, bool print) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system);

    // perform the RHF method
    res = RestrictedHartreeFock(rhfopt).run(system, ints, res, false);

    // transform the one-electron integrals to the MS basis
    ints.Tms = Transform::SingleSpin(ints.T, res.rhf.C);
    ints.Vms = Transform::SingleSpin(ints.V, res.rhf.C);

    // transform the Coulomb integral to the MS basis
    ints.Jms = Transform::CoulombSpin(ints.J, res.rhf.C);

    // run the RCI and return the total energy
    return {run(system, ints, res, print), ints};
}

Result RestrictedConfigurationInteraction::run(const System& system, const Integrals& ints, Result res, bool print) const {
    // start the timer and define the determinant vector
    Timer::Timepoint start = Timer::Now(); std::vector<Determinant> dets;

    // generate the determinants
    if (opt.excitations.size() == 0) dets = Determinant(ints.S.rows(), system.nocc(), system.nocc()).full();
    else throw std::runtime_error("RCI EXCITATIONS NOT IMPLEMENTED YET");

    // print the elapsed time
    if (print) std::cout << "GENERATED " << dets.size() << " DETERMINANTS: " << std::flush;
    if (print) std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl;

    // print the number of determinants
    if (print) std::cout << "\nFILLING RCI HAMILTONIAN: " << std::flush;

    // create the CI Hamiltonian and reset the timer
    res.rci.F = Matrix<>(dets.size(), dets.size()); start = Timer::Now();

    // fill the CI Hamiltonian
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < res.rci.F.rows(); i++) {
        for (int j = 0; j < i + 1; j++) {
            res.rci.F(i, j) = dets.at(i).hamilton(dets.at(j), ints.Tms + ints.Vms, ints.Jms); res.rci.F(j, i) = res.rci.F(i, j);
        }
    }

    // print the matrix creation time and reset the timer
    if (print) {std::cout << Timer::Format(Timer::Elapsed(start)) << "\nFINDING THE EIGENVALUES: " << std::flush;} start = Timer::Now();

    // create the eigenvalue solver
    Eigen::SelfAdjointEigenSolver<Matrix<>> solver(res.rci.F);

    // extract the eigenvectors and eigenvalues
    res.rci.C = solver.eigenvectors(); res.rci.eps = solver.eigenvalues().array() + system.repulsion();

    // print the eigenproblem time
    if (print) std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl << std::endl;

    // assign energy and return results
    res.Etot = res.rci.eps(0); return res;
}

#include "restrictedconfigurationinteraction.h"

Result RestrictedConfigurationInteraction::run(const System& system, Result res, bool print) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system);

    // perform the RHF method
    res = RestrictedHartreeFock(rhfopt).run(system, ints, res, false);

    // transform the Coulomb integral to the MO basis
    ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C);

    // run the RCI and return the total energy
    return run(system, ints, res, print);
}

Result RestrictedConfigurationInteraction::run(const System& system, const Integrals& ints, Result res, bool) const {
    return res;
}

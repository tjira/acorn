#include "restrictedmollerplesset.h"

Result RestrictedMollerPlesset::run(const System& system, Result res, bool print) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system);

    // perform the RHF method
    res = RestrictedHartreeFock(rhfopt).run(system, ints, res, false);

    // transform the Coulomb integral to the MO basis
    ints.Jmo = Transform::Coulomb(ints.J, res.rhf.C);

    // run the RMP and return the total energy
    return run(system, ints, res, print);
}

Result RestrictedMollerPlesset::run(const System& system, const Integrals& ints, Result res, bool) const {
    // initialize the correlation energy
    res.rmp.Ecorr = 0;

    // calculate the correlation energy
    for (int i = 0; i < system.nocc(); i++) {
        for (int j = 0; j < system.nocc(); j++) {
            for (int a = system.nocc(); a < ints.Jmo.dimension(0); a++) {
                for (int b = system.nocc(); b < ints.Jmo.dimension(1); b++) {
                    res.rmp.Ecorr += ints.Jmo(i, a, j, b) * (2 * ints.Jmo(i, a, j, b) - ints.Jmo(i, b, j, a)) / (res.rhf.eps(i) + res.rhf.eps(j) - res.rhf.eps(a) - res.rhf.eps(b));
                }
            }
        }
    }

    // assign the total energy and return the struct
    res.Etot = res.rhf.E + res.rmp.Ecorr + system.repulsion(); return res;
}

#include "method.cpp"
template class Method<RestrictedMollerPlesset>;

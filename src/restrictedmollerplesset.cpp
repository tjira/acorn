#include "restrictedmollerplesset.h"

Method::Result RestrictedMollerPlesset::run(const System& system, const Integrals& ints, Result res, bool) const {
    // transform the integrals to the MO basis
    Tensor<> Jmo = Transform::Coulomb(ints.J, res.rhf.C);

    // calculate the correlation energy
    for (int i = 0; i < system.nocc(); i++) {
        for (int j = 0; j < system.nocc(); j++) {
            for (int a = system.nocc(); a < Jmo.dimension(0); a++) {
                for (int b = system.nocc(); b < Jmo.dimension(1); b++) {
                    res.rmp.Ecorr += Jmo(i, a, j, b) * (2 * Jmo(i, a, j, b) - Jmo(i, b, j, a)) / (res.rhf.eps(i) + res.rhf.eps(j) - res.rhf.eps(a) - res.rhf.eps(b));
                }
            }
        }
    }

    // assign the total energy and return the struct
    res.Etot = res.rhf.E + res.rmp.Ecorr + system.repulsion(); return res;

}

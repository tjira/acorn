#include "restrictedmollerplesset.h"

std::tuple<Result, Integrals> RestrictedMollerPlesset::run(const System& system, Result res, bool print) const {
    // define the integral struct
    Integrals ints;

    // calculate all the atomic integrals
    ints.S = Integral::Overlap(system), ints.T = Integral::Kinetic(system);
    ints.V = Integral::Nuclear(system), ints.J = Integral::Coulomb(system);

    // perform the RHF method
    res = RestrictedHartreeFock(rhfopt).run(system, ints, res, false);

    // transform the Coulomb integral to the MS basis
    ints.Jms = Transform::CoulombSpin(ints.J, res.rhf.C);

    // run the RMP and return the total energy
    return {run(system, ints, res, print), ints};
}

Result RestrictedMollerPlesset::run(const System& system, const Integrals& ints, Result res, bool) const {
    // initialize the correlation energy and antisymmetrized Coulomb integrals
    res.rmp.Ecorr = 0; Tensor<4> Jmsa = ints.Jms - ints.Jms.shuffle(Eigen::array<int, 4>{0, 3, 2, 1});

    // calculate the number of virtual spatial orbitals
    int nvirt = ints.Jms.dimension(0) / 2 - system.nocc();

    // define the Coulomb integrals slice indices
    Eigen::array<Eigen::Index, 4> vvvvis = {2 * system.nocc(), 2 * system.nocc(), 2 * system.nocc(), 2 * system.nocc()}, vvvvie = {2 * nvirt, 2 * nvirt, 2 * nvirt, 2 * nvirt};
    Eigen::array<Eigen::Index, 4> ovovis = {0, 2 * system.nocc(), 0, 2 * system.nocc()}, ovovie = {2 * system.nocc(), 2 * nvirt, 2 * system.nocc(), 2 * nvirt};
    Eigen::array<Eigen::Index, 4> oovvis = {0, 0, 2 * system.nocc(), 2 * system.nocc()}, oovvie = {2 * system.nocc(), 2 * system.nocc(), 2 * nvirt, 2 * nvirt};
    Eigen::array<Eigen::Index, 4> oooois = {0, 0, 0, 0}, ooooie = {2 * system.nocc(), 2 * system.nocc(), 2 * system.nocc(), 2 * system.nocc()};

    // calculate the Coulomb integral slices and define the amplitude tensor
    Tensor<4> oooo = Jmsa.slice(oooois, ooooie), oovv = Jmsa.slice(oovvis, oovvie), ovov = Jmsa.slice(ovovis, ovovie), vvvv = Jmsa.slice(vvvvis, vvvvie); Tensor<4> t = ovov;

    // divide the amplitudes by the orbital energy differences
    for (int i = 0; i < 2 * system.nocc(); i++) {
        for (int j = 0; j < 2 * system.nocc(); j++) {
            for (int a = 0; a < 2 * nvirt; a++) {
                for (int b = 0; b < 2 * nvirt; b++) {
                    t(i, a, j, b) /= res.rhf.eps(system.nocc() + a / 2) + res.rhf.eps(system.nocc() + b / 2) - res.rhf.eps(i / 2) - res.rhf.eps(j / 2);
                }
            }
        }
    }

    // calculate MP2 correlation
    if (opt.order > 1) {
        // define the contraction axes
        Eigen::IndexPair<int> a1(0, 0), a2(1, 1), a3(2, 2), a4(3, 3);

        // add the MP2 correlation energy
        res.rmp.Ecorr -= 0.25 * Tensor<0>(ovov.contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{a1, a2, a3, a4}))(0);
    }

    // calculate MP3 correlation
    if (opt.order > 2) {
        // define the contraction axes
        Eigen::IndexPair<int> a11(1, 0), a12(3, 2), a21(0, 0), a22(1, 2), a23(2, 1), a24(3, 3);
        Eigen::IndexPair<int> b11(0, 1), b12(2, 3), b21(0, 1), b22(1, 3), b23(2, 0), b24(3, 2);
        Eigen::IndexPair<int> c11(0, 1), c12(3, 2), c21(0, 1), c22(1, 2), c23(2, 0), c24(3, 3);

        // add the MP3 correlation energy
        res.rmp.Ecorr += 0.125 * Tensor<0>(t.contract(vvvv, Eigen::array<Eigen::IndexPair<int>, 2>{a11, a12}).contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{a21, a22, a23, a24}))(0);
        res.rmp.Ecorr += 0.125 * Tensor<0>(t.contract(oooo, Eigen::array<Eigen::IndexPair<int>, 2>{b11, b12}).contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{b21, b22, b23, b24}))(0);
        res.rmp.Ecorr -= Tensor<0>(t.contract(oovv, Eigen::array<Eigen::IndexPair<int>, 2>{c11, c12}).contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{c21, c22, c23, c24}))(0);
    }

    // assign the total energy and return the struct
    res.Etot = res.rhf.E + res.rmp.Ecorr + system.repulsion(); res.Eexc = Vector<>(1); res.Eexc << res.Etot; return res;
}

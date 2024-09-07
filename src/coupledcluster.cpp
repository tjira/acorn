#include "coupledcluster.h"

torch::Tensor CoupledCluster::ccd_amplitudes(const System& system, const torch::Tensor& J_MS_AP, const torch::Tensor& E_MS_D, const torch::Tensor& T2) {
    // define the slices
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // define the linear terms
    torch::Tensor lccd1 = 0.5 * torch::einsum("abcd,cdij->abij", {J_MS_AP.index({v, v, v, v}), T2});
    torch::Tensor lccd2 = 0.5 * torch::einsum("klij,abkl->abij", {J_MS_AP.index({o, o, o, o}), T2});
    torch::Tensor lccd3 =       torch::einsum("akic,bcjk->abij", {J_MS_AP.index({v, o, o, v}), T2});

    // combine the linear terms
    lccd3 = lccd3 - lccd3.swapaxes(0, 1) - lccd3.swapaxes(2, 3) + lccd3.swapaxes(0, 1).swapaxes(2, 3);

    // define the quadratic terms
    torch::Tensor ccd1 = -0.50 * torch::einsum("klcd,acij,bdkl->abij", {J_MS_AP.index({o, o, v, v}), T2, T2});
    torch::Tensor ccd2 = -0.50 * torch::einsum("klcd,abik,cdjl->abij", {J_MS_AP.index({o, o, v, v}), T2, T2});
    torch::Tensor ccd3 =  0.25 * torch::einsum("klcd,cdij,abkl->abij", {J_MS_AP.index({o, o, v, v}), T2, T2});
    torch::Tensor ccd4 =         torch::einsum("klcd,acik,bdjl->abij", {J_MS_AP.index({o, o, v, v}), T2, T2});

    // combine the quadratic terms
    ccd1 = ccd1 - ccd1.swapaxes(0, 1), ccd2 = ccd2 - ccd2.swapaxes(2, 3), ccd4 = ccd4 - ccd4.swapaxes(2, 3);

    // return the CCD amplitudes
    return E_MS_D * (J_MS_AP.index({v, v, o, o}) + lccd1 + lccd2 + lccd3 + ccd1 + ccd2 + ccd3 + ccd4);
}

std::tuple<torch::Tensor, torch::Tensor> CoupledCluster::ccsd_amplitudes(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP, const torch::Tensor& E_MS_S, const torch::Tensor& E_MS_D, const torch::Tensor& T1, const torch::Tensor& T2) {
    // define the slices
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // define the tau and tilde tau tensors
    torch::Tensor ttau = T2 + 0.5 * torch::einsum("ai,bj->abij", {T1, T1}) - 0.5 * torch::einsum("ai,bj->abij", {T1, T1}).swapaxes(2, 3);
    torch::Tensor tau  = T2 +       torch::einsum("ai,bj->abij", {T1, T1}) -       torch::einsum("ai,bj->abij", {T1, T1}).swapaxes(2, 3);

    // define the F intermediates
    torch::Tensor Fae = (1 - torch::eye(system.virtual_spinorbitals(), torch::kDouble)) * F_MS.index({v, v}) - 0.5 * torch::einsum("me,am->ae",     {F_MS.index({o, v}),          T1  })
                                                                                                             +       torch::einsum("mafe,fm->ae",   {J_MS_AP.index({o, v, v, v}), T1  })
                                                                                                             - 0.5 * torch::einsum("mnef,afmn->ae", {J_MS_AP.index({o, o, v, v}), ttau});
    torch::Tensor Fmi = (1 - torch::eye(system.electrons(),            torch::kDouble)) * F_MS.index({o, o}) + 0.5 * torch::einsum("me,ei->mi",     {F_MS.index({o, v}),          T1  })
                                                                                                             +       torch::einsum("mnie,en->mi",   {J_MS_AP.index({o, o, o, v}), T1  })
                                                                                                             + 0.5 * torch::einsum("mnef,efin->mi", {J_MS_AP.index({o, o, v, v}), ttau});
    torch::Tensor Fme =                                                                   F_MS.index({o, v}) +       torch::einsum("mnef,fn->me",   {J_MS_AP.index({o, o, v, v}), T1  });

    // define the supporting tensors for the W intermediates
    torch::Tensor P1   =            torch::einsum("ej,mnie->mnij", {T1, J_MS_AP.index({o, o, o, v})});
    torch::Tensor P2   =            torch::einsum("bm,amef->abef", {T1, J_MS_AP.index({v, o, v, v})});
    torch::Tensor T12  = 0.5 * T2 + torch::einsum("fj,bn->fbjn",   {T1, T1                         });

    // define the W intermediates
    torch::Tensor Wmnij = J_MS_AP.index({o, o, o, o}) + 0.25 * torch::einsum("efij,mnef->mnij", {tau, J_MS_AP.index({o, o, v, v})}) + P1 - P1.swapaxes(2, 3);
    torch::Tensor Wabef = J_MS_AP.index({v, v, v, v}) + 0.25 * torch::einsum("abmn,mnef->abef", {tau, J_MS_AP.index({o, o, v, v})}) - P2 + P2.swapaxes(0, 1);
    torch::Tensor Wmbej = J_MS_AP.index({o, v, v, o}) +        torch::einsum("fj,mbef->mbej",   {T1,  J_MS_AP.index({o, v, v, v})})
                                                      -        torch::einsum("bn,mnej->mbej",   {T1,  J_MS_AP.index({o, o, v, o})})
                                                      -        torch::einsum("fbjn,mnef->mbej", {T12, J_MS_AP.index({o, o, v, v})});

    // define the right hand sides of both amplitudes
    torch::Tensor RHS_T1 = F_MS.index({v, o}).clone(), RHS_T2 = J_MS_AP.index({v, v, o, o}).clone();

    // calculate the right hand sides of the single amplitudes
    RHS_T1 +=       torch::einsum("ei,ae->ai",     {T1, Fae                        });
    RHS_T1 -=       torch::einsum("am,mi->ai",     {T1, Fmi                        });
    RHS_T1 +=       torch::einsum("aeim,me->ai",   {T2, Fme                        });
    RHS_T1 -=       torch::einsum("fn,naif->ai",   {T1, J_MS_AP.index({o, v, o, v})});
    RHS_T1 -= 0.5 * torch::einsum("efim,maef->ai", {T2, J_MS_AP.index({o, v, v, v})});
    RHS_T1 -= 0.5 * torch::einsum("aemn,nmei->ai", {T2, J_MS_AP.index({o, o, v, o})});

    // calculate the right hand sides of the double amplitudes
    torch::Tensor Fmea =       torch::einsum("bm,me->be",        {T1,     Fme                        });
    torch::Tensor Fmeb =       torch::einsum("ej,me->mj",        {T1,     Fme                        });
                  P1   =       torch::einsum("aeij,be->abij",    {T2,     Fae - 0.5 * Fmea           });
                  P2   =       torch::einsum("abim,mj->abij",    {T2,     Fmi + 0.5 * Fmeb           });
    torch::Tensor P3   =       torch::einsum("aeim,mbej->abij",  {T2,     Wmbej                      });
                  P3  -=       torch::einsum("ei,am,mbej->abij", {T1, T1, J_MS_AP.index({o, v, v, o})});
    torch::Tensor P4   =       torch::einsum("ei,abej->abij",    {T1,     J_MS_AP.index({v, v, v, o})});
    torch::Tensor P5   =       torch::einsum("am,mbij->abij",    {T1,     J_MS_AP.index({o, v, o, o})});
    RHS_T2            += 0.5 * torch::einsum("abmn,mnij->abij",  {tau,    Wmnij                      });
    RHS_T2            += 0.5 * torch::einsum("efij,abef->abij",  {tau,    Wabef                      });

    // add all necessary permutations to the right hand sides of the double amplitudes
    RHS_T2 += P1.swapaxes(0, 0) - P1.swapaxes(0, 1).swapaxes(0, 0);
    RHS_T2 -= P2.swapaxes(0, 0) - P2.swapaxes(2, 3).swapaxes(0, 0);
    RHS_T2 += P3.swapaxes(0, 0) - P3.swapaxes(2, 3).swapaxes(0, 0);
    RHS_T2 -= P3.swapaxes(0, 1) - P3.swapaxes(0, 1).swapaxes(2, 3);
    RHS_T2 += P4.swapaxes(0, 0) - P4.swapaxes(2, 3).swapaxes(0, 0);
    RHS_T2 -= P5.swapaxes(0, 0) - P5.swapaxes(0, 1).swapaxes(0, 0);

    // return the amplitudes
    return {RHS_T1 * E_MS_S, RHS_T2 * E_MS_D};
}

double CoupledCluster::ccsd_energy(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP, const torch::Tensor& T1, const torch::Tensor& T2) {
    // define the CCSD energy and the slices
    double energy_ccsd = 0; auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // calculate the CCSD energy
    energy_ccsd +=        torch::einsum("ia,ai",      {F_MS.index({o, v}),          T1    }).item<double>();
    energy_ccsd += 0.25 * torch::einsum("ijab,abij",  {J_MS_AP.index({o, o, v, v}), T2    }).item<double>();
    energy_ccsd += 0.50 * torch::einsum("ijab,ai,bj", {J_MS_AP.index({o, o, v, v}), T1, T1}).item<double>();

    // return the CCSD energy
    return energy_ccsd;
}

double CoupledCluster::ccd_energy(const System& system, const torch::Tensor& J_MS_AP, const torch::Tensor& T2) {
    // define the slices
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // return the CCD energy
    return 0.25 * torch::einsum("ijab,abij", {J_MS_AP.index({o, o, v, v}), T2}).item<double>();
}

double CoupledCluster::perturbation_triple(const System& system, const torch::Tensor& J_MS_AP, const torch::Tensor& E_MS_T, const torch::Tensor& T1, const torch::Tensor& T2) {
    // define the slices
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // define the first helper tensor
    torch::Tensor P1 = torch::einsum("ai,jkbc->abcijk", {T1, J_MS_AP.index({o, o, v, v})}); torch::Tensor T3D = P1.clone();

    // add the disconnected triples
    T3D -= torch::einsum("abcijk->bacijk", P1);
    T3D -= torch::einsum("abcijk->cbaijk", P1);
    T3D -= torch::einsum("abcijk->abcjik", P1);
    T3D += torch::einsum("abcijk->bacjik", P1);
    T3D += torch::einsum("abcijk->cbajik", P1);
    T3D -= torch::einsum("abcijk->abckji", P1);
    T3D += torch::einsum("abcijk->backji", P1);
    T3D += torch::einsum("abcijk->cbakji", P1);

    // define the second helper tensor
    torch::Tensor P2 = torch::einsum("aejk,eibc->abcijk", {T2, J_MS_AP.index({v, o, v, v})}) - torch::einsum("bcim,majk->abcijk", {T2, J_MS_AP.index({o, v, o, o})}); torch::Tensor T3C = P2.clone();

    // add the connected triples
    T3C -= torch::einsum("abcijk->bacijk", P2);
    T3C -= torch::einsum("abcijk->cbaijk", P2);
    T3C -= torch::einsum("abcijk->abcjik", P2);
    T3C += torch::einsum("abcijk->bacjik", P2);
    T3C += torch::einsum("abcijk->cbajik", P2);
    T3C -= torch::einsum("abcijk->abckji", P2);
    T3C += torch::einsum("abcijk->backji", P2);
    T3C += torch::einsum("abcijk->cbakji", P2);

    // return the perturbation energy
    return (1.0 / 36.0) * torch::einsum("abcijk,abcijk", {T3C, E_MS_T * (T3C + T3D)}).item<double>();
}

std::tuple<torch::Tensor, double>  CoupledCluster::ccd_scf(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const {
    // define the slices
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // define the excitation energy tensors
    torch::Tensor E_MS_D = Transform::ExcitationEnergyFraction(F_MS, Slice(None, system.electrons()), Slice(system.electrons(), None), 2);

    // initialize the guess for amplitudes
    torch::Tensor T2 = J_MS_AP.index({v, v, o, o}) * E_MS_D;

    // initialize the energy and print the header
    double energy_ccd = 0; std::printf("\nCCD ENERGY CALCULATION\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIMER");

    // start the CCD iterations
    for (int i = 0; i < input.max_iter; i++) {

        // define the iteration timer
        Timepoint iteration_timer = Timer::Now();

        // calculate the CCD amplitudes
        T2 = ccd_amplitudes(system, J_MS_AP, E_MS_D, T2);

        // calculate the energy and save the previous value
        double energy_ccd_prev = energy_ccd; energy_ccd = ccd_energy(system, J_MS_AP, T2);

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, energy_ccd, std::abs(energy_ccd - energy_ccd_prev), Timer::Format(Timer::Elapsed(iteration_timer)).c_str());

        // break if the error is below the threshold
        if (std::abs(energy_ccd - energy_ccd_prev) < input.threshold) break;
    }

    // return the CCD energy
    return {T2, energy_ccd};
}

std::tuple<torch::Tensor, torch::Tensor, double>  CoupledCluster::ccsd_scf(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const {
    // define the slices
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None);

    // define the excitation energy tensors
    torch::Tensor E_MS_S = Transform::ExcitationEnergyFraction(F_MS, Slice(None, system.electrons()), Slice(system.electrons(), None), 1);
    torch::Tensor E_MS_D = Transform::ExcitationEnergyFraction(F_MS, Slice(None, system.electrons()), Slice(system.electrons(), None), 2);

    // initialize the guess for amplitudes
    torch::Tensor T1 = torch::zeros({system.virtual_spinorbitals(), system.electrons()}, torch::kDouble), T2 = J_MS_AP.index({v, v, o, o}) * E_MS_D;

    // initialize the energy and print the header
    double energy_ccsd = 0; std::printf("\nCCSD ENERGY CALCULATION\n%6s %20s %8s %12s\n", "ITER", "ENERGY", "|dE|", "TIMER");

    // start the CCSD iterations
    for (int i = 0; i < input.max_iter; i++) {

        // define the iteration timer
        Timepoint iteration_timer = Timer::Now();

        // calculate the CCSD amplitudes
        std::tie(T1, T2) = ccsd_amplitudes(system, F_MS, J_MS_AP, E_MS_S, E_MS_D, T1, T2);

        // calculate the energy and save the previous value
        double energy_ccsd_prev = energy_ccsd; energy_ccsd = ccsd_energy(system, F_MS, J_MS_AP, T1, T2);

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, energy_ccsd, std::abs(energy_ccsd - energy_ccsd_prev), Timer::Format(Timer::Elapsed(iteration_timer)).c_str());

        // break if the error is below the threshold
        if (std::abs(energy_ccsd - energy_ccsd_prev) < input.threshold) break;
    }

    // return the CCSD energy
    return {T1, T2, energy_ccsd};
}

std::string CoupledCluster::get_name() const {
    // initialize the name variable
    std::string name = "CC";

    // extract the method booleans
    bool ccd = input.excitation.size() == 1 && input.excitation.at(0) == 2, ccsd = input.excitation.size() == 2 && input.excitation.at(0) == 1 && input.excitation.at(1) == 2;

    // extract the perturbation booleans
    bool pt = input.perturbation.size() == 1 && input.perturbation.at(0) == 3;

    // add the method to the name
    if (ccd) name += "D"; else if (ccsd) name += "SD";

    // add the perturbation to the name
    if (pt) name += "(T)";

    // return the name
    return name;
}

double CoupledCluster::run(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const {
    // define the slices, energy and amplitude tensors
    auto o = Slice(None, system.electrons()), v = Slice(system.electrons(), None); double energy_cc = 0; torch::Tensor T1, T2;

    // extract the method booleans
    bool ccd = input.excitation.size() == 1 && input.excitation.at(0) == 2, ccsd = input.excitation.size() == 2 && input.excitation.at(0) == 1 && input.excitation.at(1) == 2;

    // extract the perturbation booleans
    bool pt = input.perturbation.size() == 1 && input.perturbation.at(0) == 3;

    // thor error if the excitation or perturbation is not supported
    if ((!ccd && !ccsd) || (input.perturbation.size() > 1 && !pt)) throw std::runtime_error("PROVIDED EXCITATION/PERTRUBATION COMBO IS NOT SUPPORTED");

    // calculate the energy
    if (ccsd) std::tie(T1, T2, energy_cc) = ccsd_scf(system, F_MS, J_MS_AP);
    if (ccd ) std::tie(    T2, energy_cc) = ccd_scf (system, F_MS, J_MS_AP);

    // make the T1 amplitudes zero if not calculated
    if (T1.sizes().size() == 1 && T1.size(0) == 0) T1 = torch::zeros({system.virtual_spinorbitals(), system.electrons()}, torch::kDouble);

    // calculate the perturbation energy
    if (pt) energy_cc += perturbation_triple(system, J_MS_AP, Transform::ExcitationEnergyFraction(F_MS, o, v, 3), T1, T2);

    // return the energy
    return energy_cc;
}

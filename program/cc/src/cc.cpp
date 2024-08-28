#include "cc.h"

torch::Tensor Acorn::CC::LCCD::amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Emsd, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None);

    torch::Tensor lccd1 = 0.5 * torch::einsum("abcd,cdij->abij", {Jmsa.index({v, v, v, v}), T2});
    torch::Tensor lccd2 = 0.5 * torch::einsum("klij,abkl->abij", {Jmsa.index({o, o, o, o}), T2});
    torch::Tensor lccd3 =       torch::einsum("akic,bcjk->abij", {Jmsa.index({v, o, o, v}), T2});

    lccd3 = lccd3 - lccd3.swapaxes(0, 1) - lccd3.swapaxes(2, 3) + lccd3.swapaxes(0, 1).swapaxes(2, 3);

    return Emsd * (Jmsa.index({v, v, o, o}) + lccd1 + lccd2 + lccd3);
}

double Acorn::CC::LCCD::energy(const torch::Tensor& Jmsa, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None);

    torch::Tensor Elccd = 0.25 * torch::einsum("ijab,abij", {Jmsa.index({o, o, v, v}), T2});

    return Elccd.item<double>();
}

torch::Tensor Acorn::CC::CCD::amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Emsd, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None);

    torch::Tensor lccd1 = 0.5 * torch::einsum("abcd,cdij->abij", {Jmsa.index({v, v, v, v}), T2});
    torch::Tensor lccd2 = 0.5 * torch::einsum("klij,abkl->abij", {Jmsa.index({o, o, o, o}), T2});
    torch::Tensor lccd3 =       torch::einsum("akic,bcjk->abij", {Jmsa.index({v, o, o, v}), T2});

    lccd3 = lccd3 - lccd3.swapaxes(0, 1) - lccd3.swapaxes(2, 3) + lccd3.swapaxes(0, 1).swapaxes(2, 3);

    torch::Tensor ccd1 = -0.50 * torch::einsum("klcd,acij,bdkl->abij", {Jmsa.index({o, o, v, v}), T2, T2});
    torch::Tensor ccd2 = -0.50 * torch::einsum("klcd,abik,cdjl->abij", {Jmsa.index({o, o, v, v}), T2, T2});
    torch::Tensor ccd3 =  0.25 * torch::einsum("klcd,cdij,abkl->abij", {Jmsa.index({o, o, v, v}), T2, T2});
    torch::Tensor ccd4 =         torch::einsum("klcd,acik,bdjl->abij", {Jmsa.index({o, o, v, v}), T2, T2});

    ccd1 = ccd1 - ccd1.swapaxes(0, 1), ccd2 = ccd2 - ccd2.swapaxes(2, 3), ccd4 = ccd4 - ccd4.swapaxes(2, 3);

    return Emsd * (Jmsa.index({v, v, o, o}) + lccd1 + lccd2 + lccd3 + ccd1 + ccd2 + ccd3 + ccd4);
}

double Acorn::CC::CCD::energy(const torch::Tensor& Jmsa, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None);

    torch::Tensor Eccd = 0.25 * torch::einsum("ijab,abij", {Jmsa.index({o, o, v, v}), T2});

    return Eccd.item<double>();
}

std::tuple<torch::Tensor, torch::Tensor> Acorn::CC::CCSD::amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Fms, const torch::Tensor& Emss, const torch::Tensor& Emsd, const torch::Tensor& T1, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None); int nvs = Fms.sizes().at(0) - nos;

    torch::Tensor ttau = T2 + 0.5 * torch::einsum("ai,bj->abij", {T1, T1}) - 0.5 * torch::einsum("ai,bj->abij", {T1, T1}).swapaxes(2, 3);
    torch::Tensor tau  = T2 +       torch::einsum("ai,bj->abij", {T1, T1}) -       torch::einsum("ai,bj->abij", {T1, T1}).swapaxes(2, 3);

    torch::Tensor Fae = (1 - torch::eye(nvs, torch::kDouble)) * Fms.index({v, v}) - 0.5 * torch::einsum("me,am->ae",     {Fms.index({o, v}),        T1  })
                                                                                  +       torch::einsum("mafe,fm->ae",   {Jmsa.index({o, v, v, v}), T1  })
                                                                                  - 0.5 * torch::einsum("mnef,afmn->ae", {Jmsa.index({o, o, v, v}), ttau});
    torch::Tensor Fmi = (1 - torch::eye(nos, torch::kDouble)) * Fms.index({o, o}) + 0.5 * torch::einsum("me,ei->mi",     {Fms.index({o, v}),        T1  })
                                                                                  +       torch::einsum("mnie,en->mi",   {Jmsa.index({o, o, o, v}), T1  })
                                                                                  + 0.5 * torch::einsum("mnef,efin->mi", {Jmsa.index({o, o, v, v}), ttau});
    torch::Tensor Fme =                                         Fms.index({o, v}) +       torch::einsum("mnef,fn->me",   {Jmsa.index({o, o, v, v}), T1  });

    torch::Tensor Fmea =            torch::einsum("bm,me->be",   {T1, Fme});
    torch::Tensor Fmeb =            torch::einsum("ej,me->mj",   {T1, Fme});
    torch::Tensor T12  = 0.5 * T2 + torch::einsum("fj,bn->fbjn", {T1, T1 });

    torch::Tensor P1 = torch::einsum("ej,mnie->mnij", {T1, Jmsa.index({o, o, o, v})});
    torch::Tensor P2 = torch::einsum("bm,amef->abef", {T1, Jmsa.index({v, o, v, v})});

    torch::Tensor Wmnij = Jmsa.index({o, o, o, o}) + 0.25 * torch::einsum("efij,mnef->mnij", {tau, Jmsa.index({o, o, v, v})}) + P1 - P1.swapaxes(2, 3);
    torch::Tensor Wabef = Jmsa.index({v, v, v, v}) + 0.25 * torch::einsum("abmn,mnef->abef", {tau, Jmsa.index({o, o, v, v})}) - P2 + P2.swapaxes(0, 1);
    torch::Tensor Wmbej = Jmsa.index({o, v, v, o}) +        torch::einsum("fj,mbef->mbej",   {T1,  Jmsa.index({o, v, v, v})})
                                                   -        torch::einsum("bn,mnej->mbej",   {T1,  Jmsa.index({o, o, v, o})})
                                                   -        torch::einsum("fbjn,mnef->mbej", {T12, Jmsa.index({o, o, v, v})});

    torch::Tensor rhs_T1 = Fms.index({v, o}).clone(), rhs_T2 = Jmsa.index({v, v, o, o}).clone();

    rhs_T1 +=       torch::einsum("ei,ae->ai",     {T1, Fae                     });
    rhs_T1 -=       torch::einsum("am,mi->ai",     {T1, Fmi                     });
    rhs_T1 +=       torch::einsum("aeim,me->ai",   {T2, Fme                     });
    rhs_T1 -=       torch::einsum("fn,naif->ai",   {T1, Jmsa.index({o, v, o, v})});
    rhs_T1 -= 0.5 * torch::einsum("efim,maef->ai", {T2, Jmsa.index({o, v, v, v})});
    rhs_T1 -= 0.5 * torch::einsum("aemn,nmei->ai", {T2, Jmsa.index({o, o, v, o})});

                  P1  = torch::einsum("aeij,be->abij",    {T2,     Fae - 0.5 * Fmea        });
                  P2  = torch::einsum("abim,mj->abij",    {T2,     Fmi + 0.5 * Fmeb        });
    torch::Tensor P3  = torch::einsum("aeim,mbej->abij",  {T2,     Wmbej                   });
                  P3 -= torch::einsum("ei,am,mbej->abij", {T1, T1, Jmsa.index({o, v, v, o})});
    torch::Tensor P4  = torch::einsum("ei,abej->abij",    {T1,     Jmsa.index({v, v, v, o})});
    torch::Tensor P5  = torch::einsum("am,mbij->abij",    {T1,     Jmsa.index({o, v, o, o})});

    rhs_T2 += 0.5 * torch::einsum("abmn,mnij->abij", {tau, Wmnij});
    rhs_T2 += 0.5 * torch::einsum("efij,abef->abij", {tau, Wabef});
    rhs_T2 += P1.swapaxes(0, 0) - P1.swapaxes(0, 1).swapaxes(0, 0);
    rhs_T2 -= P2.swapaxes(0, 0) - P2.swapaxes(2, 3).swapaxes(0, 0);
    rhs_T2 += P3.swapaxes(0, 0) - P3.swapaxes(2, 3).swapaxes(0, 0);
    rhs_T2 -= P3.swapaxes(0, 1) - P3.swapaxes(0, 1).swapaxes(2, 3);
    rhs_T2 += P4.swapaxes(0, 0) - P4.swapaxes(2, 3).swapaxes(0, 0);
    rhs_T2 -= P5.swapaxes(0, 0) - P5.swapaxes(0, 1).swapaxes(0, 0);

    return {rhs_T1 * Emss, rhs_T2 * Emsd};
}

double Acorn::CC::CCSD::energy(const torch::Tensor& Jmsa, const torch::Tensor& Fms, const torch::Tensor& T1, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None);

    torch::Tensor Eccsd = torch::einsum("ia,ai", {Fms.index({o, v}), T1}) + 0.25 * torch::einsum("ijab,abij", {Jmsa.index({o, o, v, v}), T2}) + 0.5 * torch::einsum("ijab,ai,bj", {Jmsa.index({o, o, v, v}), T1, T1});

    return Eccsd.item<double>();
}

double Acorn::CC::CCSD::perturbationTriple(const torch::Tensor& Jmsa, const torch::Tensor& Emst, const torch::Tensor& T1, const torch::Tensor& T2, int nos) {
    auto o = Slice(None, nos), v = Slice(nos, None);

    torch::Tensor P1 = torch::einsum("ai,jkbc->abcijk", {T1, Jmsa.index({o, o, v, v})}); torch::Tensor T3D = P1.clone();

    T3D -= torch::einsum("abcijk->bacijk", P1);
    T3D -= torch::einsum("abcijk->cbaijk", P1);
    T3D -= torch::einsum("abcijk->abcjik", P1);
    T3D += torch::einsum("abcijk->bacjik", P1);
    T3D += torch::einsum("abcijk->cbajik", P1);
    T3D -= torch::einsum("abcijk->abckji", P1);
    T3D += torch::einsum("abcijk->backji", P1);
    T3D += torch::einsum("abcijk->cbakji", P1);

    torch::Tensor P2 = torch::einsum("aejk,eibc->abcijk", {T2, Jmsa.index({v, o, v, v})}) - torch::einsum("bcim,majk->abcijk", {T2, Jmsa.index({o, v, o, o})}); torch::Tensor T3C = P2.clone();

    T3C -= torch::einsum("abcijk->bacijk", P2);
    T3C -= torch::einsum("abcijk->cbaijk", P2);
    T3C -= torch::einsum("abcijk->abcjik", P2);
    T3C += torch::einsum("abcijk->bacjik", P2);
    T3C += torch::einsum("abcijk->cbajik", P2);
    T3C -= torch::einsum("abcijk->abckji", P2);
    T3C += torch::einsum("abcijk->backji", P2);
    T3C += torch::einsum("abcijk->cbakji", P2);

    return (1.0 / 36.0) * torch::einsum("abcijk,abcijk", {T3C, Emst * (T3C + T3D)}).item<double>();
}

void Acorn::CC::run(const Options& opt, std::vector<timepoint>& timers) {
    // start the timer for integral loading
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral loading
    std::cout << "SYSTEM AND INTEGRALS IN MS BASIS READING: " << std::flush;

    // load the system and integrals in MS basis from disk
    torch::Tensor Jms = torch::ReadTensor("J_MS.mat");
    torch::Tensor Vms = torch::ReadTensor("V_MS.mat");
    torch::Tensor Tms = torch::ReadTensor("T_MS.mat");
    torch::Tensor Ems = torch::ReadTensor("E_MS.mat");
    torch::Tensor Fms = torch::ReadTensor("F_MS.mat");
    torch::Tensor N   = torch::ReadTensor("N.mat"   );

    // print the time for integral loading
    std::cout << eltime(timers.at(1)) << std::endl;

    // extract the number of occupied and virtual spinorbitals
    int nos = 2 * N.index({0}).item<int>(); int nvs = Ems.sizes().at(0) - nos;

    // initialize the energy variables and orbital slices
    double E = 0, Ecc = 0, Eccp = 0; auto o = Slice(None, nos), v = Slice(nos, None);

    // initialize the antisymmetrized Coulomb integrals in Physicists' notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); torch::Tensor Hms = Vms + Tms;

    // create orbital energy tensors
    torch::Tensor Emst = 1 / (Ems.index({o}).reshape({-1}) + Ems.index({o}).reshape({-1, 1}) + Ems.index({o}).reshape({-1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1, 1, 1}));
    torch::Tensor Emsd = 1 / (Ems.index({o}).reshape({-1}) + Ems.index({o}).reshape({-1, 1}) - Ems.index({v}).reshape({-1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1}));
    torch::Tensor Emss = 1 / (Ems.index({o}).reshape({-1}) - Ems.index({v}).reshape({-1, 1}));

    // initialize single and double excitation amplitudes
    torch::Tensor T1 = torch::zeros({nvs, nos}, torch::dtype(torch::kDouble)), T2 = Jmsa.index({v, v, o, o}) * Emsd;

    // calculate the HF energy
    for (int i = 0; i < nos; i++) {
        E += Hms.index({i, i}).item<double>(); for (int j = 0; j < nos; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item<double>();
    }

    // print the required number of contributions
    std::printf("\nCC ENERGY CALCULATION\n%13s %17s %12s\n", "CONTR", "VALUE", "TIME");

    // update amplitudes to self consistency
    for (int i = 0; i < opt.iters; i++) {

        // start the timer
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // update the amplitudes
        if (std::find(opt.exc.begin(), opt.exc.end(), 1) != opt.exc.end() && std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 2 && !opt.linear) {
            std::tie(T1, T2) = Acorn::CC::CCSD::amplitude(Jmsa, Fms, Emss, Emsd, T1, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::CCSD::energy(Jmsa, Fms, T1, T2, nos);
        } else if (std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 1 && !opt.linear) {
            T2 = Acorn::CC::CCD::amplitude(Jmsa, Emsd, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::CCD::energy(Jmsa, T2, nos);
        } else if (std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 1 && opt.linear) {
            T2 = Acorn::CC::LCCD::amplitude(Jmsa, Emsd, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::LCCD::energy(Jmsa, T2, nos);
        } else {
            throw std::runtime_error("NO VALID EXCITATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
        }

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, Ecc, std::abs(Ecc - Eccp), eltime(timers.at(1)).c_str());

        // finish if covergence reached
        if (std::abs(Ecc - Eccp) < opt.thresh) {std::cout << std::endl; break;}
        else if (i == opt.iters - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE CC AMPLITUDE SCF REACHED");
    }

    // calculate the perturbation corrections
    if (std::find(opt.exc.begin(), opt.exc.end(), 1) != opt.exc.end() && std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 2 && !opt.linear) {
        if (std::find(opt.pts.begin(), opt.pts.end(), 3) != opt.pts.end() && opt.pts.size() == 1) {
            timers.at(1) = std::chrono::high_resolution_clock().now(); std::cout << "THIRD ORDER EXCITATIONS CORRECTION: " << std::flush; E += Acorn::CC::CCSD::perturbationTriple(Jmsa, Emst, T1, T2, nos); std::cout << eltime(timers.at(1)) << std::endl << std::endl;
        } else if (opt.pts.size() > 1) {
            throw std::runtime_error("NO VALID PERTURBATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
        }
    } else if (opt.pts.size() > 0) {
        throw std::runtime_error("NO VALID EXCITATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
    }

    // print the final energy
    std::printf("FINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + Ecc + N.index({1}).item<double>(), eltime(timers.at(0)).c_str());
}

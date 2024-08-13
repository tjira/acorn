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

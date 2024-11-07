---
title: Code Solutions
parent: Supplementary Material
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Code Solutions<!--\label{sec:code_solutions}-->

This section provides the solutions to all of the coding exercises provided in the text. The solutions are written in Python and use the NumPy library for numerical operations. The code snippets are self-contained and can be run in any Python environment. The solutions are organized by the exercise they correspond to and are presented in the same order as in the text. For convenience, the full code solution file [resmet.py](/acorn/python/resmet.py) can be saved here.

## Hartree--Fock Method<!--\label{sec:hf_code_solution}-->

<!--{id=code:hf_solution caption="Hartree--Fock method exercise code solution."}-->
```python
# define energies, number of occupied and virtual orbitals and the number of basis functions
E_HF, E_HF_P, VNN, nocc, nvirt, nbf = 0, 1, 0, sum(atoms) // 2, S.shape[0] - sum(atoms) // 2, S.shape[0]

# define some matrices and tensors
K, F, D, C, eps = J.transpose(0, 3, 2, 1), np.zeros_like(S), np.zeros_like(S), np.zeros_like(S), np.array(nbf * [0])

# calculate the X matrix which is the inverse of the square root of the overlap matrix
SEP = np.linalg.eigh(S); X = SEP[1] @ np.diag(1 / np.sqrt(SEP[0])) @ SEP[1].T

# the scf loop
while abs(E_HF - E_HF_P) > args.threshold:

    # build the Fock matrix
    F = H + np.einsum("ijkl,ij->kl", J - 0.5 * K, D, optimize=True)

    # solve the Fock equations
    eps, C = np.linalg.eigh(X @ F @ X); C = X @ C

    # build the density from coefficients
    D = 2 * np.einsum("ij,kj->ik", C[:, :nocc], C[:, :nocc])

    # save the previous energy and calculate the current electron energy
    E_HF_P, E_HF = E_HF, 0.5 * np.einsum("ij,ij", D, H + F, optimize=True)

# calculate nuclear-nuclear repulsion
for i, j in it.product(range(len(atoms)), range(len(atoms))):
    VNN += 0.5 * atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :]) if i != j else 0

# print the results
print("    RHF ENERGY: {:.8f}".format(E_HF + VNN))
```

## Integral Transform<!--\label{sec:int_code_solution}-->

<!--{id=code:int_solution caption="Integral transform exercise code solution."}-->
```python
# define the occ and virt spinorbital slices shorthand
o, v = slice(0, 2 * nocc), slice(2 * nocc, 2 * nbf)

# define the tiling matrix for the Jmsa coefficients and energy placeholders
P = np.array([np.eye(nbf)[:, i // 2] for i in range(2 * nbf)]).T

# define the spin masks
M = np.repeat([1 - np.arange(2 * nbf, dtype=int) % 2], nbf, axis=0)
N = np.repeat([    np.arange(2 * nbf, dtype=int) % 2], nbf, axis=0)

# tile the coefficient matrix, apply the spin mask and tile the orbital energies
Cms, epsms = np.block([[C @ P], [C @ P]]) * np.block([[M], [N]]), np.repeat(eps, 2)

# transform the core Hamiltonian and Fock matrix to the molecular spinorbital basis
Hms = np.einsum("ip,ij,jq->pq", Cms, np.kron(np.eye(2), H), Cms, optimize=True)
Fms = np.einsum("ip,ij,jq->pq", Cms, np.kron(np.eye(2), F), Cms, optimize=True)

# transform the coulomb integrals to the MS basis in chemist's notation
Jms = np.einsum("ip,jq,ijkl,kr,ls->pqrs", Cms, Cms, np.kron(np.eye(2), np.kron(np.eye(2), J).T), Cms, Cms, optimize=True);

# antisymmetrized two-electron integrals in physicist's notation
Jmsa = (Jms - Jms.swapaxes(1, 3)).transpose(0, 2, 1, 3)

# create the tensors of reciprocal differences of orbital energies in MS basis used in post-HF methods
Emss, Emsd = (epsms[o].reshape(-1) - epsms[v].reshape(-1, 1)), (epsms[o].reshape(-1) + epsms[o].reshape(-1, 1) - epsms[v].reshape(-1, 1, 1) - epsms[v].reshape(-1, 1, 1, 1))
```

## 2nd and 3rd Order Perturbative Corrections<!--\label{sec:mp_code_solution}-->

<!--{id=code:mp_solution caption="Møller--Plesset Perturbation Theory of 2nd Order and Møller--Plesset Perturbation Theory of 3rd Order exercise code solution."}-->
```python
# energy containers
E_MP2, E_MP3 = 0, 0

# calculate the MP2 correlation energy
if args.mp2 or args.mp3:
    E_MP2 += 0.25 * np.einsum("abij,ijab,abij", Jmsa[v, v, o, o], Jmsa[o, o, v, v], Emsd**-1, optimize=True)
    print("    MP2 ENERGY: {:.8f}".format(E_HF + E_MP2 + VNN))

# calculate the MP3 correlation energy
if args.mp3:
    E_MP3 += 0.125 * np.einsum("abij,cdab,ijcd,abij,cdij", Jmsa[v, v, o, o], Jmsa[v, v, v, v], Jmsa[o, o, v, v], Emsd**-1, Emsd**-1, optimize=True)
    E_MP3 += 0.125 * np.einsum("abij,ijkl,klab,abij,abkl", Jmsa[v, v, o, o], Jmsa[o, o, o, o], Jmsa[o, o, v, v], Emsd**-1, Emsd**-1, optimize=True)
    E_MP3 += 1.000 * np.einsum("abij,cjkb,ikac,abij,acik", Jmsa[v, v, o, o], Jmsa[v, o, o, v], Jmsa[o, o, v, v], Emsd**-1, Emsd**-1, optimize=True)
    print("    MP3 ENERGY: {:.8f}".format(E_HF + E_MP2 + E_MP3 + VNN))
```

## Full Configuration Interaction<!--\label{sec:ci_code_solution}-->

<!--{id=code:ci_solution caption="Configuration Interaction exercise code solution."}-->
```python
# generate the determiants
dets = [np.array(det) for det in it.combinations(range(2 * nbf), 2 * nocc)]

# define the CI Hamiltonian
Hci = np.zeros([len(dets), len(dets)])

# define the Slater-Condon rules, "so" is an array of unique and common spinorbitals [unique, common]
slater0 = lambda so: sum([Hms[m, m] for m in so]) + sum([0.5 * Jmsa[m, n, m, n] for m, n in it.product(so, so)])
slater1 = lambda so: Hms[so[0], so[1]] + sum([Jmsa[so[0], m, so[1], m] for m in so[2:]])
slater2 = lambda so: Jmsa[so[0], so[1], so[2], so[3]]

# filling of the CI Hamiltonian
for i in range(0, Hci.shape[0]):
    for j in range(i, Hci.shape[1]):

        # aligned determinant and the sign
        aligned, sign = dets[j].copy(), 1

        # align the determinant "j" to "i" and calculate the sign
        for k in (k for k in range(len(aligned)) if aligned[k] != dets[i][k]):
            while len(l := np.where(dets[i] == aligned[k])[0]) and l[0] != k:
                aligned[[k, l[0]]] = aligned[[l[0], k]]; sign *= -1

        # find the unique and common spinorbitals
        so = np.block([
            np.array([aligned[k] for k in range(len(aligned)) if aligned[k] not in dets[i]]),
            np.array([dets[i][k] for k in range(len(dets[j])) if dets[i][k] not in aligned]),
            np.array([aligned[k] for k in range(len(aligned)) if aligned[k] in dets[i]])
        ]).astype(int)

        # apply the Slater-Condon rules and multiply by the sign
        if ((aligned - dets[i]) != 0).sum() == 0: Hci[i, j] = slater0(so) * sign
        if ((aligned - dets[i]) != 0).sum() == 1: Hci[i, j] = slater1(so) * sign
        if ((aligned - dets[i]) != 0).sum() == 2: Hci[i, j] = slater2(so) * sign

        # fill the lower triangle
        Hci[j, i] = Hci[i, j]

# solve the eigensystem and assign energy
eci, Cci = np.linalg.eigh(Hci); E_FCI = eci[0] - E_HF

# print the results
print("    FCI ENERGY: {:.8f}".format(E_HF + E_FCI + VNN))
```

## Coupled Cluster Singles and Doubles<!--\label{sec:cc_code_solution}-->

<!--{id=code:cc_solution caption="Coupled Cluster Doubles and Coupled Cluster Singles and Doubles method exercise code solution."}-->
```python
# energy containers for all the CC methods
E_CCD, E_CCD_P, E_CCSD, E_CCSD_P = 0, 1, 0, 1

# initialize the first guess for the t-amplitudes as zeros
t1, t2 = np.zeros((2 * nvirt, 2 * nocc)), np.zeros((2 * nvirt, 2 * nvirt, 2 * nocc, 2 * nocc))

# CCD loop
if args.ccd:
    while abs(E_CCD - E_CCD_P) > args.threshold:

        # collect all the distinct LCCD terms
        lccd1 = 0.5 * np.einsum("abcd,cdij->abij", Jmsa[v, v, v, v], t2, optimize=True)
        lccd2 = 0.5 * np.einsum("klij,abkl->abij", Jmsa[o, o, o, o], t2, optimize=True)
        lccd3 =       np.einsum("akic,bcjk->abij", Jmsa[v, o, o, v], t2, optimize=True)

        # apply the permuation operator and add it to the corresponding LCCD terms
        lccd3 = lccd3 + lccd3.transpose(1, 0, 3, 2) - lccd3.transpose(1, 0, 2, 3) - lccd3.transpose(0, 1, 3, 2)

        # collect all the distinct first CCD terms
        ccd1 = -0.50 * np.einsum("klcd,acij,bdkl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)
        ccd2 = -0.50 * np.einsum("klcd,abik,cdjl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)
        ccd3 =  0.25 * np.einsum("klcd,cdij,abkl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)
        ccd4 =         np.einsum("klcd,acik,bdjl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)

        # apply the permuation operator and add it to the corresponding CCD terms
        ccd1, ccd2, ccd4 = ccd1 - ccd1.transpose(1, 0, 2, 3), ccd2 - ccd2.transpose(0, 1, 3, 2), ccd4 - ccd4.transpose(0, 1, 3, 2)

        # update the t-amplitudes
        t2 = (Jmsa[v, v, o, o] + lccd1 + lccd2 + lccd3 + ccd1 + ccd2 + ccd3 + ccd4) / Emsd

        # evaluate the energy
        E_CCD_P, E_CCD = E_CCD, 0.25 * np.einsum("ijab,abij", Jmsa[o, o, v, v], t2, optimize=True)

    # print the CCD energy
    print("    CCD ENERGY: {:.8f}".format(E_HF + E_CCD + VNN))

# CCSD loop
if args.ccsd:
    while abs(E_CCSD - E_CCSD_P) > args.threshold:

        # calculate the effective two-particle excitation operators
        ttau = t2 + 0.5 * np.einsum("ai,bj->abij", t1, t1, optimize=True) - 0.5 * np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(2, 3)
        tau  = t2 +       np.einsum("ai,bj->abij", t1, t1, optimize=True) -       np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(2, 3)

        # calculate the 2D two-particle intermediates
        Fae = (1 - np.eye(2 * nvirt)) * Fms[v, v] - 0.5 * np.einsum("me,am->ae",     Fms[o, v],        t1,   optimize=True)                                                   +       np.einsum("mafe,fm->ae",   Jmsa[o, v, v, v], t1,   optimize=True)                                                   - 0.5 * np.einsum("mnef,afmn->ae", Jmsa[o, o, v, v], ttau, optimize=True)
        Fmi = (1 - np.eye(2 * nocc )) * Fms[o, o] + 0.5 * np.einsum("me,ei->mi",     Fms[o, v],        t1,   optimize=True)                                                   +       np.einsum("mnie,en->mi",   Jmsa[o, o, o, v], t1,   optimize=True)                                                   + 0.5 * np.einsum("mnef,efin->mi", Jmsa[o, o, v, v], ttau, optimize=True)
        Fme =                           Fms[o, v] +       np.einsum("mnef,fn->me",   Jmsa[o, o, v, v], t1,   optimize=True)

        # define some complementary variables used in the following expressions
        Fmea =            np.einsum("bm,me->be",   t1, Fme, optimize=True)
        Fmeb =            np.einsum("ej,me->mj",   t1, Fme, optimize=True)
        t12  = 0.5 * t2 + np.einsum("fj,bn->fbjn", t1, t1,  optimize=True)

        # define the permutation arguments for all terms the W intermediates
        P1 = np.einsum("ej,mnie->mnij", t1, Jmsa[o, o, o, v], optimize=True)
        P2 = np.einsum("bm,amef->abef", t1, Jmsa[v, o, v, v], optimize=True)

        # calculate the 4D two-particle intermediates
        Wmnij = Jmsa[o, o, o, o] + 0.25 * np.einsum("efij,mnef->mnij", tau, Jmsa[o, o, v, v], optimize=True) + P1 - P1.swapaxes(2, 3)
        Wabef = Jmsa[v, v, v, v] + 0.25 * np.einsum("abmn,mnef->abef", tau, Jmsa[o, o, v, v], optimize=True) - P2 + P2.swapaxes(0, 1)
        Wmbej = Jmsa[o, v, v, o] +        np.einsum("fj,mbef->mbej",   t1,  Jmsa[o, v, v, v], optimize=True)                                  -        np.einsum("bn,mnej->mbej",   t1,  Jmsa[o, o, v, o], optimize=True)                                  -        np.einsum("fbjn,mnef->mbej", t12, Jmsa[o, o, v, v], optimize=True)

        # define the right hand side of the T1 and T2 amplitude equations
        rhs_T1, rhs_T2 = Fms[v, o].copy(), Jmsa[v, v, o, o].copy()

        # calculate the right hand side of the CCSD equation for T1
        rhs_T1 +=       np.einsum("ei,ae->ai",     t1, Fae,              optimize=True)
        rhs_T1 -=       np.einsum("am,mi->ai",     t1, Fmi,              optimize=True)
        rhs_T1 +=       np.einsum("aeim,me->ai",   t2, Fme,              optimize=True)
        rhs_T1 -=       np.einsum("fn,naif->ai",   t1, Jmsa[o, v, o, v], optimize=True)
        rhs_T1 -= 0.5 * np.einsum("efim,maef->ai", t2, Jmsa[o, v, v, v], optimize=True)
        rhs_T1 -= 0.5 * np.einsum("aemn,nmei->ai", t2, Jmsa[o, o, v, o], optimize=True)

        # define the permutation arguments for all terms in the equation for T2
        P1  = np.einsum("aeij,be->abij",    t2,     Fae - 0.5 * Fmea, optimize=True)
        P2  = np.einsum("abim,mj->abij",    t2,     Fmi + 0.5 * Fmeb, optimize=True)
        P3  = np.einsum("aeim,mbej->abij",  t2,     Wmbej,            optimize=True)
        P3 -= np.einsum("ei,am,mbej->abij", t1, t1, Jmsa[o, v, v, o], optimize=True)
        P4  = np.einsum("ei,abej->abij",    t1,     Jmsa[v, v, v, o], optimize=True)
        P5  = np.einsum("am,mbij->abij",    t1,     Jmsa[o, v, o, o], optimize=True)

        # calculate the right hand side of the CCSD equation for T2
        rhs_T2 += 0.5 * np.einsum("abmn,mnij->abij", tau, Wmnij, optimize=True)
        rhs_T2 += 0.5 * np.einsum("efij,abef->abij", tau, Wabef, optimize=True)
        rhs_T2 += P1.transpose(0, 1, 2, 3) - P1.transpose(1, 0, 2, 3)
        rhs_T2 -= P2.transpose(0, 1, 2, 3) - P2.transpose(0, 1, 3, 2)
        rhs_T2 += P3.transpose(0, 1, 2, 3) - P3.transpose(0, 1, 3, 2)
        rhs_T2 -= P3.transpose(1, 0, 2, 3) - P3.transpose(1, 0, 3, 2)
        rhs_T2 += P4.transpose(0, 1, 2, 3) - P4.transpose(0, 1, 3, 2)
        rhs_T2 -= P5.transpose(0, 1, 2, 3) - P5.transpose(1, 0, 2, 3)

        # Update T1 and T2 amplitudes
        t1, t2 = rhs_T1 / Emss, rhs_T2 / Emsd

        # evaluate the energy
        E_CCSD_P, E_CCSD = E_CCSD, np.einsum("ia,ai", Fms[o, v], t1) + 0.25 * np.einsum("ijab,abij", Jmsa[o, o, v, v], t2) + 0.5 * np.einsum("ijab,ai,bj", Jmsa[o, o, v, v], t1, t1)

    # print the CCSD energy
    print("   CCSD ENERGY: {:.8f}".format(E_HF + E_CCSD + VNN))
```

{:.cite}

---
title: Code Solutions
parent: Supplementary Material
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Code Solutions<!--\label{sec:code_solutions}-->

This section provides the solutions to all of the coding exercises provided in the text. The solutions are written in Python and use the NumPy library for numerical operations. The code snippets are self-contained and can be run in any Python environment. The solutions are organized by the exercise they correspond to and are presented in the same order as in the text. For convenience, the full code solution files [resmet.py](/acorn/python/resmet.py), [neqdet.py](/acorn/python/neqdet.py) and [shcdet.py](/acorn/python/shcdet.py) can be saved by clicking the links.

## Restricted Electronic Structure Methods<!--\label{sec:resmet_code}-->

<!--{id=code:resmet caption="Restricted electronic structure methods code."}-->
```python
#!/usr/bin/env python

import argparse as ap, itertools as it, numpy as np

np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%20.14f" % x))

ATOM = {
    "H" :  1,                                                                                                                                                                 "He":  2,
    "Li":  3, "Be":  4,                                                                                                     "B" :  5, "C" :  6, "N" :  7, "O" :  8, "F" :  9, "Ne": 10,
    "Na": 11, "Mg": 12,                                                                                                     "Al": 13, "Si": 14, "P" : 15, "S" : 16, "Cl": 17, "Ar": 18,
    "K" : 19, "Ca": 20, "Sc": 21, "Ti": 22, "V" : 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y" : 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I" : 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Hf": 72, "Ta": 73, "W" : 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88, "Ac": 89, "Rf":104, "Db":105, "Sg":106, "Bh":107, "Hs":108, "Mt":109, "Ds":110, "Rg":111, "Cn":112, "Nh":113, "Fl":114, "Mc":115, "Lv":116, "Ts":117, "Og":118,

                                  "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
                                  "Th": 90, "Pa": 91, "U" : 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm":100, "Md":101, "No":102, "Lr":103,
}

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =========================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="resmet.py", description="Restricted Electronic Structure Methods Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-m", "--molecule", help="Molecule file in the .xyz format. (default: %(default)s)", type=str, default="molecule.xyz")
    parser.add_argument("-t", "--threshold", help="Convergence threshold for SCF loop. (default: %(default)s)", type=float, default=1e-12)
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)

    # add integral options
    parser.add_argument("--int", help="Filenames of the integral files. (default: %(default)s)", nargs=4, type=str, default=["T_AO.mat", "V_AO.mat", "S_AO.mat", "J_AO.mat"])

    # optional flags
    parser.add_argument("--diis", help="Enable the DIIS convergence accelerator.", action=ap.BooleanOptionalAction)

    # method switches
    parser.add_argument("--fci", help="Perform the full configuration interaction calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--ccd", help="Perform the doubles coupled clusters calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--ccsd", help="Perform the singles/doubles coupled clusters calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--mp2", help="Perform the second order Moller-Plesset calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--mp3", help="Perform the third order Moller-Plesset calculation.", action=ap.BooleanOptionalAction)

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # forward declaration of integrals and energies in MS basis (just to avoid NameError in linting)
    Hms, Fms, Jms, Jmsa, Emss, Emsd = np.zeros(2 * [0]), np.zeros(2 * [0]), np.zeros(4 * [0]), np.zeros(4 * [0]), np.array([[]]), np.array([[[[]]]])

    # OBTAIN THE MOLECULE AND ATOMIC INTEGRALS ===========================================

    # get the atomic numbers and coordinates of all atoms
    atoms  = np.array([ATOM[line.split()[0]] for line in open(args.molecule).readlines()[2:]], dtype=int  )
    coords = np.array([line.split()[1:]      for line in open(args.molecule).readlines()[2:]], dtype=float)

    # convert coordinates to bohrs and forward declare orbital slices and atom count
    coords *= 1.8897261254578281; o, v = slice(0, 0), slice(0, 0); natoms = len(atoms)

    # load the integrals from the files
    H, S = np.loadtxt(args.int[0], skiprows=1) + np.loadtxt(args.int[1], skiprows=1), np.loadtxt(args.int[2], skiprows=1); J = np.loadtxt(args.int[3], skiprows=1).reshape(4 * [S.shape[1]])

    # HARTREE-FOCK METHOD ================================================================

    # energies, iteration counter, number of basis functions and number of occupied/virtual orbitals
    E_HF, E_HF_P, VNN, iter, nbf, nocc = 0, 1, 0, 0, S.shape[0], sum(atoms) // 2; nvirt = nbf - nocc

    # exchange integrals and the guess density matrix
    K, D = J.transpose(0, 3, 2, 1), np.zeros((nbf, nbf))

    # Fock matrix, coefficient matrix and orbital energies initialized to zero
    F, C, eps = np.zeros((nbf, nbf)), np.zeros((nbf, nbf)), np.zeros((nbf))

    # DIIS containers
    DIIS_F, DIIS_E = [], []

    # the X matrix which is the inverse of the square root of the overlap matrix
    SEP = np.linalg.eigh(S); X = SEP[1] @ np.diag(1 / np.sqrt(SEP[0])) @ SEP[1].T

    # the scf loop
    while abs(E_HF - E_HF_P) > args.threshold:

        # build the Fock matrix and increment the iteration counter
        F = H + np.einsum("ijkl,ij->kl", J - 0.5 * K, D, optimize=True); iter += 1

        # DIIS extrapolation
        if args.diis and iter > 1:

            # append the DIIS matrices
            DIIS_F.append(F); DIIS_E.append(S @ D @ F - F @ D @ S);

            # truncate the DIIS subspace
            if len(DIIS_F) > 5: DIIS_F.pop(0), DIIS_E.pop(0)

            # build the DIIS system
            A = np.ones ((len(DIIS_F) + 1, len(DIIS_F) + 1)); A[-1, -1] = 0
            b = np.zeros((len(DIIS_F) + 1                 )); b[-1]     = 1

            # fill the DIIS matrix
            for i, j in it.product(range(len(DIIS_F)), range(len(DIIS_F))):
                A[i, j] = A[j, i] = np.einsum("ij,ij", DIIS_E[i], DIIS_E[j])

            # solve the DIIS equations and extrapolate the Fock matrix
            c = np.linalg.solve(A, b); F = np.einsum("i,ijk->jk", c[:-1], DIIS_F)

        # solve the Fock equations
        eps, C = np.linalg.eigh(X @ F @ X); C = X @ C

        # build the density from coefficients
        D = 2 * np.einsum("ij,kj->ik", C[:, :nocc], C[:, :nocc])

        # save the previous energy and calculate the current electron energy
        E_HF_P, E_HF = E_HF, 0.5 * np.einsum("ij,ij", D, H + F, optimize=True)

    # calculate nuclear-nuclear repulsion
    for i, j in ((i, j) for i, j in it.product(range(natoms), range(natoms)) if i != j):
        VNN += 0.5 * atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :])

    # print the results
    print("    RHF ENERGY: {:.8f} ({} ITERATIONS)".format(E_HF + VNN, iter))

    # INTEGRAL TRANSFORMS FOR POST-HARTREE-FOCK METHODS ==================================
    if args.mp2 or args.mp3 or args.ccd or args.ccsd or args.fci:

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

        # transform the coulomb integrals to the MS basis in chemists' notation
        Jms = np.einsum("ip,jq,ijkl,kr,ls->pqrs",
            Cms, Cms, np.kron(np.eye(2), np.kron(np.eye(2), J).T), Cms, Cms, optimize=True
        );

        # antisymmetrized two-electron integrals in physicists' notation
        Jmsa = (Jms - Jms.swapaxes(1, 3)).transpose(0, 2, 1, 3)

        # tensor epsilon_i^a
        Emss = epsms[o] - epsms[v, None]

        # tensor epsilon_ij^ab
        Emsd = epsms[o] + epsms[o, None] - epsms[v, None, None] - epsms[v, None, None, None]

    # MOLLER-PLESSET PERTURBATION THEORY =================================================
    if args.mp2 or args.mp3:

        # energy containers
        E_MP2, E_MP3 = 0, 0

        # calculate the MP2 correlation energy
        if args.mp2 or args.mp3:
            E_MP2 += 0.25 * np.einsum("abij,ijab,abij",
                Jmsa[v, v, o, o], Jmsa[o, o, v, v], Emsd**-1, optimize=True
            )
            print("    MP2 ENERGY: {:.8f}".format(E_HF + E_MP2 + VNN))

        # calculate the MP3 correlation energy
        if args.mp3:
            E_MP3 += 0.125 * np.einsum("abij,cdab,ijcd,abij,cdij",
                Jmsa[v, v, o, o], Jmsa[v, v, v, v], Jmsa[o, o, v, v], Emsd**-1, Emsd**-1,
                optimize=True
            )
            E_MP3 += 0.125 * np.einsum("abij,ijkl,klab,abij,abkl",
                Jmsa[v, v, o, o], Jmsa[o, o, o, o], Jmsa[o, o, v, v], Emsd**-1, Emsd**-1,
                optimize=True
            )
            E_MP3 += 1.000 * np.einsum("abij,cjkb,ikac,abij,acik",
                Jmsa[v, v, o, o], Jmsa[v, o, o, v], Jmsa[o, o, v, v], Emsd**-1, Emsd**-1,
                optimize=True
            )
            print("    MP3 ENERGY: {:.8f}".format(E_HF + E_MP2 + E_MP3 + VNN))

    # COUPLED CLUSTER METHOD =============================================================
    if args.ccd or args.ccsd:

        # energy containers for all the CC methods
        E_CCD, E_CCD_P, E_CCSD, E_CCSD_P = 0, 1, 0, 1

        # initialize the first guess for the t-amplitudes as zeros
        t1, t2 = np.zeros((2 * nvirt, 2 * nocc)), np.zeros(2 * [2 * nvirt] + 2 * [2 * nocc])

        # CCD loop
        if args.ccd:
            while abs(E_CCD - E_CCD_P) > args.threshold:

                # collect all the distinct LCCD terms
                lccd1 = 0.5 * np.einsum("abcd,cdij->abij", Jmsa[v, v, v, v], t2, optimize=True)
                lccd2 = 0.5 * np.einsum("klij,abkl->abij", Jmsa[o, o, o, o], t2, optimize=True)
                lccd3 = 1.0 * np.einsum("akic,bcjk->abij", Jmsa[v, o, o, v], t2, optimize=True)

                # apply the permuation operator and add it to the corresponding LCCD terms
                lccd3 += lccd3.transpose(1, 0, 3, 2) - lccd3.swapaxes(0, 1) - lccd3.swapaxes(2, 3)

                # collect all the remaining CCD terms
                ccd1 = -0.50 * np.einsum("klcd,acij,bdkl->abij",
                    Jmsa[o, o, v, v], t2, t2, optimize=True
                )
                ccd2 = -0.50 * np.einsum("klcd,abik,cdjl->abij",
                    Jmsa[o, o, v, v], t2, t2, optimize=True
                )
                ccd3 = +0.25 * np.einsum("klcd,cdij,abkl->abij",
                    Jmsa[o, o, v, v], t2, t2, optimize=True
                )
                ccd4 = +1.00 * np.einsum("klcd,acik,bdjl->abij",
                    Jmsa[o, o, v, v], t2, t2, optimize=True
                )

                # permutation operators
                ccd1 -= ccd1.swapaxes(0, 1);
                ccd2 -= ccd2.swapaxes(2, 3);
                ccd4 -= ccd4.swapaxes(2, 3)

                # update the t-amplitudes
                t2 = (Jmsa[v, v, o, o] + lccd1 + lccd2 + lccd3 + ccd1 + ccd2 + ccd3 + ccd4) / Emsd

                # evaluate the energy
                E_CCD_P, E_CCD = E_CCD, 0.25 * np.einsum("ijab,abij", Jmsa[o, o, v, v], t2)

            # print the CCD energy
            print("    CCD ENERGY: {:.8f}".format(E_HF + E_CCD + VNN))

        # CCSD loop
        if args.ccsd:
            while abs(E_CCSD - E_CCSD_P) > args.threshold:

                # define taus
                tau, ttau = t2.copy(), t2.copy()

                # add contributions to the tilde tau
                ttau += 0.5 * np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(0, 0)
                ttau -= 0.5 * np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(2, 3)

                # add the contributions to tau
                tau += np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(0, 0)
                tau -= np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(2, 3)

                # define the deltas for Fae and Fmi
                dae, dmi = np.eye(2 * nvirt), np.eye(2 * nocc)

                # define Fae, Fmi and Fme
                Fae, Fmi, Fme = (1 - dae) * Fms[v, v], (1 - dmi) * Fms[o, o], Fms[o, v].copy()

                # add the contributions to Fae
                Fae -= 0.5 * np.einsum("me,am->ae",     Fms[o, v],        t1,   optimize=True)
                Fae += 1.0 * np.einsum("mafe,fm->ae",   Jmsa[o, v, v, v], t1,   optimize=True)
                Fae -= 0.5 * np.einsum("mnef,afmn->ae", Jmsa[o, o, v, v], ttau, optimize=True)

                # add the contributions to Fmi
                Fmi += 0.5 * np.einsum("me,ei->mi",     Fms[o, v],        t1,   optimize=True)
                Fmi += 1.0 * np.einsum("mnie,en->mi",   Jmsa[o, o, o, v], t1,   optimize=True)
                Fmi += 0.5 * np.einsum("mnef,efin->mi", Jmsa[o, o, v, v], ttau, optimize=True)

                # add the contributions to Fme
                Fme += np.einsum("mnef,fn->me", Jmsa[o, o, v, v], t1, optimize=True)

                # define Wmnij, Wabef and Wmbej
                Wmnij = Jmsa[o, o, o, o].copy()
                Wabef = Jmsa[v, v, v, v].copy()
                Wmbej = Jmsa[o, v, v, o].copy()

                # define some complementary variables used in the Wmbej intermediate
                t12  = 0.5 * t2 + np.einsum("fj,bn->fbjn", t1, t1,  optimize=True)

                # add contributions to Wmnij
                Wmnij += 0.25 * np.einsum("efij,mnef->mnij", tau, Jmsa[o, o, v, v], optimize=True)
                Wabef += 0.25 * np.einsum("abmn,mnef->abef", tau, Jmsa[o, o, v, v], optimize=True)
                Wmbej += 1.00 * np.einsum("fj,mbef->mbej",   t1,  Jmsa[o, v, v, v], optimize=True) 
                Wmbej -= 1.00 * np.einsum("bn,mnej->mbej",   t1,  Jmsa[o, o, v, o], optimize=True) 
                Wmbej -= 1.00 * np.einsum("fbjn,mnef->mbej", t12, Jmsa[o, o, v, v], optimize=True) 

                # define the permutation arguments for Wmnij and Wabef and add them
                PWmnij = np.einsum("ej,mnie->mnij", t1, Jmsa[o, o, o, v], optimize=True)
                PWabef = np.einsum("bm,amef->abef", t1, Jmsa[v, o, v, v], optimize=True)

                # add the permutations to Wmnij and Wabef
                Wmnij += PWmnij - PWmnij.swapaxes(2, 3)
                Wabef += PWabef.swapaxes(0, 1) - PWabef

                # define the right hand side of the T1 and T2 amplitude equations
                rhs_T1, rhs_T2 = Fms[v, o].copy(), Jmsa[v, v, o, o].copy()

                # calculate the right hand side of the CCSD equation for T1
                rhs_T1 += 1.0 * np.einsum("ei,ae->ai",     t1, Fae,              optimize=True)
                rhs_T1 -= 1.0 * np.einsum("am,mi->ai",     t1, Fmi,              optimize=True)
                rhs_T1 += 1.0 * np.einsum("aeim,me->ai",   t2, Fme,              optimize=True)
                rhs_T1 -= 1.0 * np.einsum("fn,naif->ai",   t1, Jmsa[o, v, o, v], optimize=True)
                rhs_T1 -= 0.5 * np.einsum("efim,maef->ai", t2, Jmsa[o, v, v, v], optimize=True)
                rhs_T1 -= 0.5 * np.einsum("aemn,nmei->ai", t2, Jmsa[o, o, v, o], optimize=True)

                # contracted F matrices that used in the T2 equations
                Faet = np.einsum("bm,me->be", t1, Fme, optimize=True)
                Fmet = np.einsum("ej,me->mj", t1, Fme, optimize=True)

                # define the permutation arguments for all terms in the equation for T2
                P1  = np.einsum("aeij,be->abij",    t2,     Fae - 0.5 * Faet, optimize=True)
                P2  = np.einsum("abim,mj->abij",    t2,     Fmi + 0.5 * Fmet, optimize=True)
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

                # Update T1 and T2 amplitudes and save the previous iteration
                t1, t2 = rhs_T1 / Emss, rhs_T2 / Emsd; E_CCSD_P = E_CCSD

                # evaluate the energy
                E_CCSD  = 1.00 * np.einsum("ia,ai",      Fms[o, v],        t1    )
                E_CCSD += 0.25 * np.einsum("ijab,abij",  Jmsa[o, o, v, v], t2    )
                E_CCSD += 0.50 * np.einsum("ijab,ai,bj", Jmsa[o, o, v, v], t1, t1)

            # print the CCSD energy
            print("   CCSD ENERGY: {:.8f}".format(E_HF + E_CCSD + VNN))

    # CONFIGURATION INTERACTION ==========================================================
    if args.fci:

        # generate the determiants
        dets = [np.array(det) for det in it.combinations(range(2 * nbf), 2 * nocc)]

        # define the CI Hamiltonian
        Hci = np.zeros([len(dets), len(dets)])

        # define the Slater-Condon rules, "so" is an array of unique and common spinorbitals
        slater0 = lambda so: (
            sum(np.diag(Hms)[so]) + sum([0.5 * Jmsa[m, n, m, n] for m, n in it.product(so, so)])
        )
        slater1 = lambda so: (
            Hms[so[0], so[1]] + sum([Jmsa[so[0], m, so[1], m] for m in so[2:]])
        )
        slater2 = lambda so: (
            Jmsa[so[0], so[1], so[2], so[3]]
        )

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
                so = np.block(list(map(lambda l: np.array(l), [
                    [aligned[k] for k in range(len(aligned)) if aligned[k] not in dets[i]],
                    [dets[i][k] for k in range(len(dets[j])) if dets[i][k] not in aligned],
                    [aligned[k] for k in range(len(aligned)) if aligned[k] in dets[i]]
                ]))).astype(int)

                # apply the Slater-Condon rules and multiply by the sign
                if ((aligned - dets[i]) != 0).sum() == 0: Hci[i, j] = Hci[j, i] = slater0(so) * sign
                if ((aligned - dets[i]) != 0).sum() == 1: Hci[i, j] = Hci[j, i] = slater1(so) * sign
                if ((aligned - dets[i]) != 0).sum() == 2: Hci[i, j] = Hci[j, i] = slater2(so) * sign

        # solve the eigensystem and assign energy
        eci, Cci = np.linalg.eigh(Hci); E_FCI = eci[0] - E_HF

        # print the results
        print("    FCI ENERGY: {:.8f}".format(E_HF + E_FCI + VNN))
```

## Nonadiabatic Exact Quantum Dynamics<!--\label{sec:neqdet_code}-->

<!--{id=code:neqdet caption="Nonadiabatic exact quantum dynamics."}-->
```python
#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, scipy as sp, scipy.linalg; np.random.seed(0)

# EXAMPLES
# ./neqdet.py -g "[np.exp(-(r1-1)**2-(r2-1)**2)]" -v "[[0.5*(r1**2+r2**2)]]"
# ./neqdet.py -d 0 -f 0.01 -p 8192 -l 32 -m 2000 -t 10 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" -g "[0,np.exp(-(r1+10)**2+10j*r1)]" --align --adiabatic

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =========================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="resmet.py", description="Numerically Exact Quantum Dynamics Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-a", "--align", help="Align the wavefunction plot to the potential.", action=ap.BooleanOptionalAction)
    parser.add_argument("-c", "--imaginary", help="Perform imaginary-time propagation and find specified number of orthogonal states.", type=int, default=0)
    parser.add_argument("-d", "--damp", help="Gaussian damping parameter for autocorrelation function. (default: %(default)s)", type=float, default=0.003)
    parser.add_argument("-f", "--factor", help="Factor to scale the wavefunction before plotting. (default: %(default)s)", type=float, default=10)
    parser.add_argument("-g", "--guess", help="Initial guess for the wavefunction. (default: %(default)s)", type=str, default="[np.exp(-(r1-1)**2)]")
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
    parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=500)
    parser.add_argument("-l", "--limit", help="Distance limit of the wavefunction grid. (default: %(default)s)", type=float, default=8)
    parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=float, default=1)
    parser.add_argument("-n", "--ntraj", help="Number of Bohmian trajectories. (default: %(default)s)", type=int, default=100)
    parser.add_argument("-p", "--points", help="Number of points in the wavefunction grid. (default: %(default)s)", type=int, default=128)
    parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-u", "--adiabatic", help="Transform the results to the adiabatic basis. (default: %(default)s)", action=ap.BooleanOptionalAction)
    parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")

    # parse the arguments
    args = parser.parse_args()

    # replace the variables in the string expressions
    for i in range(1, 10): args.guess     =     args.guess.replace(f"r{i}", f"r[:, {i - 1}]")
    for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # PREPARE VARIABLES FOR THE DYNAMICS =================================================

    # extract the number of dimensions without using regex
    ndim = max([int(e[0]) for e in args.potential.split("r[:, ")[1:]]) + 1

    # define the function to evaluate potential
    potf = lambda: np.array(eval(args.potential)).transpose(2, 0, 1);

    # define the function to evaluate the guess wavefunction on a specified grid
    psif = lambda: np.array(list(map(lambda x: x * np.ones(args.points**ndim), eval(args.guess))), dtype=complex).T

    # REAL AND IMAGINARY QUANTUM DYNAMICS ================================================

    # calculate the space step
    dr = 2 * args.limit / (args.points - 1)

    # create the grid in real and fourier space
    r = np.stack(np.meshgrid(*[np.linspace(-args.limit, args.limit, args.points)] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)
    k = np.stack(np.meshgrid(*[2 * np.pi * np.fft.fftfreq(args.points, dr)      ] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)

    # create the potential, wavefunction and extract wavefunction dimensions
    V, psi = potf(), psif(); shape = ndim * [args.points] + [psi.shape[1]]

    # normalize the guess wavefunction
    psi /= np.sqrt(dr * np.einsum("ij,ij->", psi.conj(), psi))

    # create the containers for the observables and wavefunctions
    position, momentum, ekin, epot = [], [], [], []; wfnopt, wfn = [], [psi]

    # iterate over the propagations
    for i in range(args.imaginary if args.imaginary else 1):

        # print the propagation header
        print() if i else None; print("PROPAGATION OF STATE %d " % (i))

        # create the initial points for Bohmian trajectories
        trajs = np.concatenate((r[np.random.choice(psi.shape[0], size=args.ntraj, p=(np.abs(psi)**2 * dr).sum(axis=1))][:, None, :], np.zeros((args.ntraj, args.iterations, ndim))), axis=1)

        # initialize the initial wavefunction from the wfn container and clear all containers
        psi = wfn[0]; wfn.clear(), position.clear(), momentum.clear(), ekin.clear(), epot.clear()

        # calculate the propagators for each point in the grid
        K = np.array([sp.linalg.expm(-0.5 * (1 if args.imaginary else 1j) * args.timestep * np.sum(k[i, :] ** 2) / args.mass * np.eye(psi.shape[1])) for i in range(r.shape[0])])
        R = np.array([sp.linalg.expm(-0.5 * (1 if args.imaginary else 1j) * args.timestep                                    * V[i]                ) for i in range(r.shape[0])])

        # print the propagation header
        print("%6s %12s %12s %12s" % ("ITER", "EKIN", "EPOT", "ETOT", ))

        # propagate the wavefunction
        for j in range(args.iterations + 1):

            # propagate in real space 
            if (j): psi = np.einsum("ijk,ik->ij", R, psi)

            # fourier transform the wavefunction
            if (j): psi = np.fft.fftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape)

            # propagate in momentum space
            if (j): psi = np.einsum("ijk,ik->ij", K, psi)

            # inverse fourier transform the wavefunction
            if (j): psi = np.fft.ifftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape)

            # propagate in real space
            if (j): psi = np.einsum("ijk,ik->ij", R, psi)

            # orthogonalize the wavefunction
            for i in range(len(wfnopt)): psi -= np.sum(wfnopt[i].conj() * psi) * wfnopt[i] * dr

            # normalize the wavefunction
            if args.imaginary: psi /= np.sqrt(dr * np.einsum("ij,ij->", psi.conj(), psi))

            # append the potential energy and the wavefunction
            epot.append(np.einsum("ij,ijk,ik->", psi.conj(), V, psi).real * dr); wfn.append(psi.copy())

            # create a n dimensional copy of the wavefunction, its fourier transform and a container for Bohmian trajectory velocity
            psid, psik, v = psi.reshape(shape), np.fft.fftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape), np.zeros(shape[:-1] + [ndim])

            # append the kinetic energy
            ekin.append((psi.conj() * np.fft.ifftn((psik * (0.5 * np.sum(k**2, axis=1) / args.mass)[:, None]).reshape(shape), axes=range(ndim)).reshape(psi.shape)).real.sum() * dr)

            # append the position
            position.append(np.sum(r * np.sum(np.abs(psi)**2, axis=1, keepdims=True), axis=0) * dr)

            # append the momentum
            momentum.append(np.array([np.sum(psi.conj() * np.fft.ifftn((1j * (k[:, dim:dim + 1]) * psik).reshape(shape), axes=range(ndim)).reshape(psi.shape)).imag * dr for dim in range(ndim)]))

            # calculate the velocity of the Bohmian trajectories
            if (j): v[..., :] = np.array([(np.conjugate(psid) * np.gradient(psid, dr, axis=dim)).sum(axis=-1).imag / ((np.abs(psid)**2).sum(axis=-1) + 1e-14) / args.mass for dim in range(ndim)]).T

            # propagate the Bohmian trajectories
            if (j): trajs[:, j, :] = trajs[:, j - 1, :] + sp.interpolate.interpn(points=ndim * [np.unique(r[:, 0])], values=v, xi=trajs[:, j - 1, :]) * args.timestep

            # print the iteration info
            if j % 100 == 0: print("%6d %12.6f %12.6f %12.6f" % (j, ekin[-1], epot[-1], ekin[-1] + epot[-1]))

        # append the optimized wavefunction to the container
        if args.imaginary: wfnopt.append(psi.copy())

    # ADIABATIC TRANSFORM AND SPECTRUM ===================================================

    # calculate the adiabatic eigenstates
    U = [np.linalg.eigh(V[i])[1] for i in range(V.shape[0])]

    # adiabatize the potential and wavefunctions
    V   = np.einsum("ijk,ikl,ilm->ijm", U, V,   U) if args.adiabatic else np.array(V  )
    wfn = np.einsum("ikl,jik->jil",     U, wfn   ) if args.adiabatic else np.array(wfn)

    # calculate the density matrices and acf
    density = np.einsum("jia,jib->jab", wfn,           wfn.conj()).real * dr
    acf     = np.einsum("ij,tij->t",    wfn[0].conj(), wfn       )      * dr

    # symmetrize the acf and apply the damping function
    acf = np.concatenate((np.flip(acf)[:-1], np.array(acf).conj())) * np.exp(-args.damp * (np.arange(-args.iterations, args.iterations + 1) * args.timestep)**2)

    # calculate the spectrum of the zero-padded acf and the corresponding energies
    spectrum = np.abs(np.fft.fft(np.pad(acf, 2 * [10 * len(acf)], mode="constant")))**2; omega = 2 * np.pi * np.fft.fftfreq(len(spectrum), args.timestep)

    # PRINT AND PLOT THE RESULTS =========================================================

    # print the final population
    print("
FINAL POPULATION: %s" % np.diag(density[-1]))

    # scale the wavefunction and add the potential to the wavefunction
    wfn = args.factor * np.array(wfn); wfn += (1 + 1j) * np.einsum("ijj,k->kij", V, np.ones(args.iterations + 1)) if args.align else 0

    # create the stationary subplots and scale the wavefunction
    fig, axs = plt.subplots(2, 3, figsize=(12, 6))

    # plot the energies
    axs[0, 0].plot(np.arange(args.iterations + 1) * args.timestep, ekin, label="Kinetic Energy"  )
    axs[0, 0].plot(np.arange(args.iterations + 1) * args.timestep, epot, label="Potential Energy")

    # plot the population
    axs[0, 1].plot(np.arange(args.iterations + 1) * args.timestep, [np.diag(rho) for rho in density], label=[f"S$_{i}$" for i in range(density.shape[1])])

    # print the acf
    axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).real, label="Re(ACF)")
    axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).imag, label="Im(ACF)")

    # plot the spectrum
    axs[1, 1].plot(omega[np.argsort(omega)], spectrum[np.argsort(omega)] / np.max(spectrum))

    # plot the position and momentum of the Bohmian trajectories
    [axs[0, 2].plot(np.arange(args.iterations + 1) * args.timestep, trajs[i, :, 0],                                                 alpha=0.05, color="tab:blue") for i in range(args.ntraj)]
    [axs[1, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.gradient(trajs[i, :, 0], args.timestep, axis=0) * args.mass, alpha=0.05, color="tab:blue") for i in range(args.ntraj)]

    # plot the numerically exact expectation values of position and momentum
    axs[0, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.array(position)[:, 0], color="tab:orange", label="$<\Psi|\hat{r_x}|\Psi>$")
    axs[1, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.array(momentum)[:, 0], color="tab:orange", label="$<\Psi|\hat{p_x}|\Psi>$")

    # plot the mean of the Bohmian trajectories
    axs[0, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.mean(trajs[:, :, 0], axis=0),                                                 "--", color="black", label="$<r_x>$")
    axs[1, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.mean(np.gradient(trajs[:, :, 0], args.timestep, axis=1), axis=0) * args.mass, "--", color="black", label="$<p_x>$")

    # set the labels for the stationary plot
    axs[0, 0].set_xlabel("Time (a.u.)"    ); axs[0, 0].set_ylabel("Energy (a.u.)"           )
    axs[0, 1].set_xlabel("Time (a.u.)"    ); axs[0, 1].set_ylabel("Population"              )
    axs[0, 2].set_xlabel("Time (a.u.)"    ); axs[0, 2].set_ylabel("Position (a.u.)"         )
    axs[1, 0].set_xlabel("Time (a.u.)"    ); axs[1, 0].set_ylabel("Autocorrelation Function")
    axs[1, 1].set_xlabel("Energy (a.u.)"  ); axs[1, 1].set_ylabel("Normalized Intensity"    )
    axs[1, 2].set_xlabel("Time (a.u.)"    ); axs[1, 2].set_ylabel("Momentum (a.u.)"         )

    # set the domain for the spectrum plot, the end will be as last element from the end less than some value
    axs[1, 1].set_xlim(0, omega[np.argsort(omega)][np.where(spectrum[np.argsort(omega)] / np.max(spectrum) > 1e-6)][-1])

    # enable legends
    axs[0, 0].legend(); axs[0, 1].legend(); axs[1, 0].legend(); axs[0, 2].legend(); axs[1, 2].legend()

    # only for 1D
    if ndim == 1:
        
        # creace the wavefunction plot
        wfig, wax = plt.subplots(1, 1)

        # plot the wavefunction
        wfnplot = np.array([[wax.plot(r, wfn[0][:, i].real, label="Re($\Psi$)")[0], wax.plot(r, wfn[0][:, i].imag, label="Im($\Psi$)")[0]] for i in range(psi.shape[1])]).flatten()

        # set the labels for the wavefunction plot and enable legend
        wax.set_xlabel("Position (a.u.)"); wax.set_ylabel("Wavefunction"); wax.legend()

        # extract the wfn min and max
        minwfn, maxwfn = min(np.array(wfn).real.min(), np.array(wfn).imag.min()), max(np.array(wfn).real.max(), np.array(wfn).imag.max())

        # set the limits for the wavefunction animation
        wax.set_ylim(minwfn - 0.1 * (maxwfn - minwfn), maxwfn + 0.1 * (maxwfn - minwfn))

        # define the update function for the wavefunction animation
        update = lambda i: [wfnplot[j].set_ydata(wfn[i][:, j // 2].real if j % 2 == 0 else wfn[i][:, j // 2].imag) for j in range(2 * psi.shape[1])]

        # make the wavefunction plot animation and set the layout
        ani = anm.FuncAnimation(wfig, update, frames=range(len(wfn)), repeat=True, interval=30); wfig.tight_layout()

    # show the plot
    fig.tight_layout(); plt.show()
```

## Surface Hopping Classical Dynamics <!--\label{sec:shcdet_code}-->

<!--{id=code:shcdet caption="Surface Hopping Classical Dynamics."}-->
```python
#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, scipy as sp, scipy.linalg

# EXAMPLES
# ./shcdet.py -p 10 1 -r -10 0.5 -s 1 -m 2000 -t 1 -i 5000 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" --adiabatic -n 1000 --lzsh

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =========================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="resmet.py", description="Surface Hopping Classical Dynamics Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
    parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=5000)
    parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=float, default=1)
    parser.add_argument("-n", "--trajectories", help="Number of classical trajectories in the simulation. (default: %(default)s)", type=int, default=1000)
    parser.add_argument("-s", "--state", help="Initial state of the trajectories. (default: %(default)s)", type=int, default=0)
    parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-p", "--momentum", help="Mean and standard deviation of momentum for each dimension. (default: %(default)s)", type=float, nargs= "+", default=[0, 1])
    parser.add_argument("-r", "--position", help="Mean and standard deviation of position for each dimension. (default: %(default)s)", type=float, nargs= "+", default=[1, 0.5])
    parser.add_argument("-u", "--adiabatic", help="Transform the results to the adiabatic basis. (default: %(default)s)", action=ap.BooleanOptionalAction)
    parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")

    # surface hopping arguments
    parser.add_argument("--lzsh", help="Enable LZSH algorithm. (default: %(default)s)", action=ap.BooleanOptionalAction)

    # parse the arguments
    args = parser.parse_args()

    # replace the variables in the string expressions
    for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # SURFACE HOPPING ALGORITHMS =========================================================

    # define the state vector and containers for the potential and transformation matrices
    s = np.zeros((args.trajectories, args.iterations + 1), dtype=int) + args.state; Vs, Us = [], []

    # define the Landau-Zener Surface Hopping algorithm
    def lzsh(i, rn):

        # loop over the trajectories and states
        for (j, k) in it.product(range(args.trajectories), range(V.shape[1])):

            # calculate the first derivatives of the new state
            dk0 = (Vs[-1][j, k, k] - Vs[-2][j, k, k]) / args.timestep
            dk1 = (Vs[-2][j, k, k] - Vs[-3][j, k, k]) / args.timestep

            # calculate the first derivatives of the current state
            ds0 = (Vs[-1][j, s[j, i], s[j, i]] - Vs[-2][j, s[j, i], s[j, i]]) / args.timestep
            ds1 = (Vs[-2][j, s[j, i], s[j, i]] - Vs[-3][j, s[j, i], s[j, i]]) / args.timestep

            # calculate the second derivatives of the new and current states
            ddk = (dk0 - dk1) / args.timestep; dds = (ds0 - ds1) / args.timestep

            # check if the trajectory is in the place for a hop
            if (k == s[j, i - 1] or (dk0 - ds0) * (dk1 - ds1) > 0 or ddk * dds > 0 or ddk - dds == 0): continue

            # calculate the hopping probability
            p = np.exp(-0.5 * np.pi * np.sqrt(sqrtarg)) if ((sqrtarg := (V[j, k, k] - V[j, s[j, i], s[j, i]])**3 / (ddk - dds)) > 0) else 0

            # hop to another state if accepted
            if (not np.isnan(p) and rn < p): s[j, i:] = k

    # PERFORM THE CLASSICAL DYNAMICS =====================================================

    # create the initial conditions for the trajectories
    r = np.column_stack([np.random.normal(loc=mu, scale=sigma, size=args.trajectories) for mu, sigma in zip(args.position[0::2], args.position[1::2])])
    v = np.column_stack([np.random.normal(loc=mu, scale=sigma, size=args.trajectories) for mu, sigma in zip(args.momentum[0::2], args.momentum[1::2])])

    # divide momentum by mass and define acceleration
    v /= args.mass; a = np.zeros_like(r); diff = 0.0001;

    # print the propagation header
    print("%6s %12s %12s %12s" % ("ITER", "EKIN", "EPOT", "ETOT", ))

    # loop over the iterations
    for i in range(args.iterations + 1):

        # save the previous values
        rp = r.copy(); vp = v.copy(); ap = a.copy()

        # loop over the dimensions and propagate
        for j in (j for j in range(r.shape[1]) if i):

            # calculate offsets
            rplus = r.copy(); rplus[:, j] += diff
            rmins = r.copy(); rmins[:, j] -= diff

            # calculate the potential energy at the offsetted positions
            Vplus = np.array(eval(args.potential.replace("r", "rplus"))).transpose(2, 0, 1)
            Vmins = np.array(eval(args.potential.replace("r", "rmins"))).transpose(2, 0, 1)

            # transform the potential energy to the adiabatic basis
            if args.adiabatic: Vplus = np.einsum("ij,jk->ijk", np.linalg.eigvalsh(Vplus), np.eye(V.shape[1]))
            if args.adiabatic: Vmins = np.einsum("ij,jk->ijk", np.linalg.eigvalsh(Vmins), np.eye(V.shape[1]))

            # calculate the acceleration
            a[:, j] = 0.5 * (np.einsum("ijj->ij", Vmins)[np.arange(args.trajectories), s[:, i]] - np.einsum("ijj->ij", Vplus)[np.arange(args.trajectories), s[:, i]]) / (diff * args.mass)

            # update the positions and velocities
            v[:, j] += 0.5 * (a[:, j] + ap[:, j]) * args.timestep; r[:, j] += (v[:, j] + 0.5 * a[:, j] * args.timestep) * args.timestep

        # calculate the diabatic potential
        V = np.array(eval(args.potential)).transpose(2, 0, 1)

        # adiabatize the potential
        E, U = np.linalg.eigh(V); V = np.einsum("ij,jk->ijk", E, np.eye(V.shape[1])); Vs.append(V); Us.append(U)

        # surface hopping algorithm
        if (args.lzsh and i > 1): lzsh(i, np.random.rand())

        # calculate the potential and kinetic energy for each trajectory
        Epot = V[np.arange(args.trajectories), s[:, i], s[:, i]]; Ekin = 0.5 * args.mass * np.sum(v**2, axis=1)

        # get the mask for trajectories that hopped to another state
        mask = (s[:, i] != s[:, i - bool(i)]) & (Ekin > (Epot - V[np.arange(args.trajectories), s[:, i - bool(i)], s[:, i - bool(i)]]))

        # rescale the velocities of the trajectories that hopped to another state
        if mask.any(): v[mask] *= np.sqrt((Ekin - Epot + V[np.arange(args.trajectories), s[:, i - 1], s[:, i - 1]]) / Ekin)[mask, None]

        # recalculate the kinetic energy if the velocities were rescaled
        if mask.any(): Ekin = 0.5 * args.mass * np.sum(v**2, axis=1)

        # print the iteration
        if i % 1000 == 0: print("%6d %12.6f %12.6f %12.6f" % (i, np.mean(Ekin), np.mean(Epot), np.mean(Ekin + Epot)))

    # PRINT AND PLOT THE RESULTS =========================================================

    # print the final populations
    print("
FINAL POPULATION: %s" % (np.bincount(s[:, -1]) / s.shape[0]))

    # create the subplots
    fig, axs = plt.subplots(1, 1, figsize=(8, 6));

    # plot the population
    axs.plot(np.arange(args.iterations + 1) * args.timestep, [np.bincount(s[:, i]) / s.shape[0] for i in range(s.shape[1])])

    # set the labels
    axs.set_xlabel("Time (a.u.)"); axs.set_ylabel("Population")

    # show the plot
    plt.tight_layout(); plt.show()
```

{:.cite}

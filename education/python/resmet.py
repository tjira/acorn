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

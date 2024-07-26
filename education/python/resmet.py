#!/usr/bin/env python

import argparse as ap, itertools as it, numpy as np

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
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

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
    parser.add_argument("--int", help="Filenames of the integral files. (default: %(default)s)", nargs=4, type=str, default=["S_AO.mat", "T_AO.mat", "V_AO.mat", "J_AO.mat"])

    # method switches
    parser.add_argument("--cisd", help="Perform the singles/doubles configuration interaction calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--fci", help="Perform the full configuration interaction calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--mp2", help="Perform the second order Moller-Plesset calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--mp3", help="Perform the third order Moller-Plesset calculation.", action=ap.BooleanOptionalAction)
    parser.add_argument("--mpn", help="Perform the n-th order Moller-Plesset calculation.", type=int, default=1)

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # forward declaration of integrals and energies in MS basis (just to avoid NameError in linting)
    Hms, Jms, Jmsa, epsms = np.zeros(2 * [0]), np.zeros(4 * [0]), np.zeros(4 * [0]), np.array([])

    # OBTAIN THE MOLECULE AND ATOMIC INTEGRALS =========================================================================================================================================================

    # get the atomic numbers and coordinates of all atoms
    atoms = np.array([ATOM[line.split()[0]] for line in open(args.molecule).readlines()[2:]], dtype=int)
    coords = np.array([line.split()[1:] for line in open(args.molecule).readlines()[2:]], dtype=float)

    # convert coordinates to bohrs
    coords *= 1.8897261254578281

    # load the integrals from the files
    S, T, V = np.loadtxt(args.int[0], skiprows=1), np.loadtxt(args.int[1], skiprows=1), np.loadtxt(args.int[2], skiprows=1); J = np.loadtxt(args.int[3], skiprows=1).reshape(4 * [S.shape[1]])

    # HARTREE-FOCK METHOD ==============================================================================================================================================================================

    # define energies, number of occupied orbitals and nbf
    E_HF, E_HF_P, VNN, nocc, nvirt, nbf = 0, 1, 0, sum(atoms) // 2, S.shape[0] - sum(atoms) // 2, S.shape[0]

    # define some matrices and tensors
    K, eps = J.transpose(0, 3, 2, 1), np.array(nbf * [0])
    H, D, C = T + V, np.zeros_like(S), np.zeros_like(S)

    # calculate the X matrix which is the inverse of the square root of the overlap matrix
    SEP = np.linalg.eigh(S); X = SEP[1] @ np.diag(1 / np.sqrt(SEP[0])) @ SEP[1].T

    while abs(E_HF - E_HF_P) > args.threshold:
        # build the Fock matrix
        F = H + np.einsum("ijkl,ij", J - 0.5 * K, D, optimize=True)

        # solve the Fock equations
        eps, C = np.linalg.eigh(X @ F @ X); C = X @ C

        # build the density from coefficients
        D = 2 * np.einsum("ij,kj", C[:, :nocc], C[:, :nocc])

        # save the previous energy and calculate the current electron energy
        E_HF_P, E_HF = E_HF, 0.5 * np.einsum("ij,ij", D, H + F, optimize=True)

    # calculate nuclear-nuclear repulsion
    for i, j in it.product(range(len(atoms)), range(len(atoms))):
        VNN += 0.5 * atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :]) if i != j else 0

    # print the results
    print("RHF ENERGY: {:.8f}".format(E_HF + VNN))

    # INTEGRAL TRANSFORMS FOR POST-HARTREE-FOCK METHODS =================================================================================================================================================
    if args.mp2 or args.mp3 or args.cisd or args.fci or args.mpn > 1:

        # define the tiling matrix for the MO coefficients and energy placeholders
        P = np.array([np.eye(nbf)[:, i // 2] for i in range(2 * nbf)]).T

        # define the spin masks
        M = np.repeat([1 - np.arange(2 * nbf, dtype=int) % 2], nbf, axis=0)
        N = np.repeat([np.arange(2 * nbf, dtype=int) % 2], nbf, axis=0)

        # tile the coefficient matrix, apply the spin mask and tile the orbital energies
        Cms, epsms = np.block([[C @ P], [C @ P]]) * np.block([[M], [N]]), np.repeat(eps, 2)

        # transform the core Hamiltonian to the molecular spinorbital basis
        Hms = np.einsum("ip,ij,jq", Cms, np.kron(np.eye(2), H), Cms, optimize=True)

        # transform the coulomb integrals to the MS basis in chemists' notation and define the antisymmetrized version
        Jms = np.einsum("ip,jq,ijkl,kr,ls", Cms, Cms, np.kron(np.eye(2), np.kron(np.eye(2), J).T), Cms, Cms, optimize=True); Jmsa = Jms - Jms.swapaxes(1, 3)

    # MOLLER-PLESSET PERTRUBATION THEORY ===============================================================================================================================================================
    if args.mp2 or args.mp3 or args.mpn > 1:

        # define the orbital slices and containers for MP2, MP3 and higher order energies
        o, v, E_MP2, E_MP3, E_MPN = slice(0, 2 * nocc), slice(2 * nocc, 2 * nbf), 0, 0, []

        # define a tensor with MP expression denominators
        epsms_vovo = 1 / (-epsms[v].reshape(-1, 1, 1, 1) + epsms[o].reshape(-1, 1, 1) - epsms[v].reshape(-1, 1) + epsms[o].reshape(-1))

        # calculate the MP2 correlation energy
        E_MP2 += 0.25 * np.einsum("aibj,iajb,aibj", Jmsa[v, o, v, o], Jmsa[o, v, o, v], epsms_vovo, optimize=True)

        # print the MP2 energy
        print("MP2 ENERGY: {:.8f}".format(E_HF + E_MP2 + VNN))

        # calculate the MP3 correlation energy
        if args.mp3 or args.mpn > 2:
            E_MP3 += 0.125 * np.einsum("aibj,ikjl,kalb,aibj,akbl", Jmsa[v, o, v, o], Jmsa[o, o, o, o], Jmsa[o, v, o, v], epsms_vovo, epsms_vovo, optimize=True)
            E_MP3 += 0.125 * np.einsum("aibj,cadb,icjd,aibj,cidj", Jmsa[v, o, v, o], Jmsa[v, v, v, v], Jmsa[o, v, o, v], epsms_vovo, epsms_vovo, optimize=True)
            E_MP3 +=         np.einsum("aibj,iack,jbkc,aibj,bjck", Jmsa[v, o, v, o], Jmsa[o, v, v, o], Jmsa[o, v, o, v], epsms_vovo, epsms_vovo, optimize=True)

        # print the MP3 energy
        if args.mp3 or args.mpn > 2: print("MP3 ENERGY: {:.8f}".format(E_HF + E_MP2 + E_MP3 + VNN))

        # function to return corresponding Coulomb integral slices
        def slices(n, contr):
            views = [Jmsa.swapaxes(1, 2)[tuple([(v if x in "abcdefgh" else o) for x in contr.split(",")[term]])] for term in range(n)]
            epsvec = [[(-epsms[v] if char in "abcdefgh" else epsms[o]) for char in contr.split(",")[i]] for i in range(n, 2 * n - 1)]
            return views + [1 / sum([v[i].reshape(-1, *[1 for _ in range(len(v) - i - 1)]) for i in range(len(v))]) for v in epsvec]
        for i in range(4, args.mpn + 1):
            exprs = [(x.split(";")[0].replace("-", ","), eval(x.split(";")[1])) for x in open(f"../../example/diagram/mbpt{i}.txt").readlines()[1:-1]]
            E_MPN.append(sum([contr[1] * np.einsum(contr[0], *slices(i, contr[0]), optimize=True) for contr in exprs]))
            print("MP{} ENERGY: {:.8f}".format(i, E_HF + E_MP2 + E_MP3 + sum(E_MPN) + VNN))

    # CONFIGUIRATION INTERACTION =======================================================================================================================================================================
    if args.cisd or args.fci:

        # define functions to generate single and double excitations of a ground configuration
        def single(ground, nbf):
            for i, j in it.product(range(len(ground)), range(len(ground), nbf)):
                yield np.array(ground[:i] + [j] + ground[i + 1:])
        def double(ground, nbf):
            for i, j in it.product(range(len(ground)), range(len(ground), nbf)):
                for k, l in it.product(range(i + 1, len(ground)), range(j + 1, nbf)):
                    yield np.array(ground[:i] + [j] + ground[i + 1:k] + [l] + ground[k + 1:])

        # define the determinant list
        dets = list()

        # generate all possible determinants for CISD
        if args.cisd:
            dets = [np.concatenate((2 * np.array(range(nocc)), 2 * np.array(range(nocc)) + 1))]
            for electron in single(list(range(nocc)), nbf):
                dets.append(np.concatenate((2 * np.array(electron), 2 * np.array(range(nocc)) + 1)))
                dets.append(np.concatenate((2 * np.array(range(nocc)), 2 * np.array(electron) + 1)))
            for alpha, beta in it.product(single(list(range(nocc)), nbf), single(list(range(nocc)), nbf)):
                dets.append(np.concatenate((2 * np.array(alpha), 2 * np.array(beta) + 1)))
            for electron in double(list(range(nocc)), nbf):
                dets.append(np.concatenate((2 * np.array(electron), 2 * np.array(range(nocc)) + 1)))
                dets.append(np.concatenate((2 * np.array(range(nocc)), 2 * np.array(electron) + 1)))

        # generate all possible determinants for FCI
        if args.fci: dets = [np.concatenate((2 * np.array(alpha), 2 * np.array(beta) + 1)) for alpha, beta in it.product(it.combinations(range(nbf), nocc), it.combinations(range(nbf), nocc))]

        # define the CI Hamiltonian
        Hci = np.zeros([len(dets), len(dets)])

        # define the Slater-Condon rules, "so" is an array of unique and common spinorbitals [unique, common]
        slater0 = lambda so: sum([Hms[m, m] for m in so]) + sum([0.5 * (Jms[m, m, n, n] - Jms[m, n, n, m]) for m, n in it.product(so, so)])
        slater1 = lambda so: Hms[so[0], so[1]] + sum([Jms[so[0], so[1], m, m] - Jms[so[0], m, m, so[1]] for m in so[2:]])
        slater2 = lambda so: Jms[so[0], so[2], so[1], so[3]] - Jms[so[0], so[3], so[1], so[2]]

        # fill the CI Hamiltonian
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
        print("{} ENERGY: {:.8f}".format("CISD" if args.cisd else "FCI", E_HF + E_FCI + VNN))

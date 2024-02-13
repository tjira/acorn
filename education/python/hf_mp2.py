#!/usr/bin/env python

import itertools as it, numpy as np

A2BOHR = 1.8897259886

ATOM = {
    "H": 1,
    "O": 8,
}

if __name__ == "__main__":
    # define some variables
    molfile, thresh = "molecule.xyz", 1e-12

    # get the atomic numbers and coordinates of all atoms
    atoms = np.array([ATOM[line.split()[0]] for line in open(molfile).readlines()[2:]], dtype=int)
    coords = np.array([line.split()[1:] for line in open(molfile).readlines()[2:]], dtype=float)

    # INTEGRAL LOADING OR CALCULATION ==================================================================================================================================================================

    # load one-electron integrals
    S, T, V = np.loadtxt("S.mat"), np.loadtxt("T.mat"), np.loadtxt("V.mat")

    # load two-electron integrals
    J = np.loadtxt("J.mat").reshape(4 * [S.shape[1]])

    # HARTREE-FOCK METHOD ==============================================================================================================================================================================

    # define energies, number of occupied orbitals and nbf
    E_HF, E_MP2, E_HF_P, VNN, nocc, nbf = 0, 0, 1, 0, sum(atoms) // 2, S.shape[0]

    # define some matrices and tensors
    H, D, C = T + V, np.zeros_like(S), np.array([])
    K, eps = J.transpose(0, 3, 2, 1), np.array([])

    # calculate the X matrix which is the inverse of the square root of the overlap matrix
    SEP = np.linalg.eigh(S); X = np.linalg.inv(SEP[1] * np.sqrt(SEP[0]) @ np.linalg.inv(SEP[1]))

    while abs(E_HF - E_HF_P) > thresh:
        # build the Fock matrix
        F = H + np.einsum("ijkl,ij", J - 0.5 * K, D, optimize=True)

        # solve the Fock equations
        eps, C = np.linalg.eigh(X @ F @ X); C = X @ C

        # build the density from coefficients
        D = 2 * np.einsum("ij,kj", C[:, :nocc], C[:, :nocc])

        # save the previous energy and calculate the current electron energy
        E_HF_P, E_HF = E_HF, 0.5 * np.einsum("ij,ij", D, H + F, optimize=True)

    # MOLLER-PLESSET PERTRUBATION THEORY OF 2ND ORDER ==================================================================================================================================================

    # transform the coulomb integral to MO basis
    Jmo = np.einsum("ip,jq,ijkl,kr,ls", C, C, J, C, C, optimize=True)

    # calculate the MP2 correlation
    for i, j, a, b in it.product(range(nocc), range(nocc), range(nocc, nbf), range(nocc, nbf)):
        E_MP2 += Jmo[i, a, j, b] * (2 * Jmo[i, a, j, b] - Jmo[i, b, j, a]) / (eps[i] + eps[j] - eps[a] - eps[b]);

    # FULL CONFIGUIRATION INTERACTION ==================================================================================================================================================================

    # create the mask of spin indices used to correctly zero out elements in tiled MO basis matrices
    spinind = np.arange(2 * nbf, dtype=int) % 2; spinmask = spinind.reshape(-1, 1) == spinind

    # transform the core Hamiltonian to the MO basis, tile it for alpha/beta spins and apply spin mask
    Hms = np.repeat(np.repeat(np.einsum("ip,ij,jq", C, H, C, optimize=True), 2, axis=0), 2, axis=1) * spinmask

    # tile the coefficient matrix to accound for different spins
    Cms = np.block([[C, np.zeros_like(C)], [np.zeros_like(C), C]])

    # transform the coulomb integral tensor to the MS basis
    Jms = np.einsum("ip,jq,ijkl,kr,ls", Cms, Cms, np.kron(np.eye(2), np.kron(np.eye(2), J).T), Cms, Cms, optimize=True)

    # NUCLEAR REPULSION AND OTHER RESULTS ==============================================================================================================================================================

    # calculate nuclear-nuclear repulsion
    for i, j in it.product(range(len(atoms)), range(len(atoms))):
        VNN += 0.5 * atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :]) / A2BOHR if i != j else 0

    # print the results
    print("HF ENERGY:", E_HF + VNN)
    print("MP2 ENERGY:", E_HF + E_MP2 + VNN)

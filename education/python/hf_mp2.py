#!/usr/bin/env python

import numpy as np

A2BOHR = 1.8897259886

ATOM = {
    "H": 1,
    "O": 8,
}

def ints():
    T, V, S = np.loadtxt("T.mat"), np.loadtxt("V.mat"), np.loadtxt("S.mat")
    J = np.loadtxt("J.mat").reshape(4 * [S.shape[1]]); return T, V, S, J

def mol(filename="molecule.xyz"):
    with open(filename) as M:
        lines = np.array([line.split() for line in M.readlines()[2:]])
    return [ATOM[S] for S in lines[:, 0]], lines[:, 1:].astype(float)

if __name__ == "__main__":
    # load the integrals and define the convergence threshold
    [T, V, S, J], [atoms, coords], thresh = ints(), mol(), 1e-12

    # define energies and number of occupied orbitals
    E_HF, E_MP2, E_HF_P, VNN, nocc = 0, 0, 1, 0, sum(atoms) // 2

    # calculate nuclear-nuclear repulsion
    for i in range(len(atoms)):
        for j in range(i):
            VNN += atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :]) / A2BOHR

    # define some matrices and tensors
    H, D = T + V, np.zeros_like(S)
    K = J.transpose(0, 3, 2, 1)

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

    # transform the coulomb integral to MO basis
    Jmo = np.einsum("ip,jq,ijkl,kr,ls", C, C, J, C, C, optimize=True)

    # calculate the MP2 correlation
    for i in range(nocc):
        for j in range(nocc):
            for a in range(nocc, len(H)):
                for b in range(nocc, len(H)):
                    E_MP2 += Jmo[i, a, j, b] * (2 * Jmo[i, a, j, b] - Jmo[i, b, j, a]) / (eps[i] + eps[j] - eps[a] - eps[b]);

    # print the results
    print("HF ENERGY:", E_HF + VNN)
    print("MP2 ENERGY:", E_HF + E_MP2 + VNN)

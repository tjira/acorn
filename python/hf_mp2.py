import numpy as np
import scipy as sp

import itertools

A2BOHR = 1.8897259886

ATOM = {
    "H": 1,
    "O": 8,
}

class IO:
    class Load:
        def Integrals():
            T, V, S = np.loadtxt("T.mat"), np.loadtxt("V.mat"), np.loadtxt("S.mat")
            J = np.loadtxt("J.mat").reshape(4 * [S.shape[0]]); return T, V, S, J
        def Molecule(M="molecule.xyz"):
            lines = np.array([line.split() for line in open(M).readlines()[2:]])
            return [ATOM[S] for S in lines[:, 0]], lines[:, 1:].astype(float)

class Determinant:
    def __init__(self, alpha, beta):
        self.alpha, self.beta = alpha, beta

    def element(self, det):
        return 1;


if __name__ == "__main__":
    # load the integrals and molecule and define the convergence threshold
    [T, V, S, J], [atoms, coords], thresh = IO.Load.Integrals(), IO.Load.Molecule(), 1e-12

    # define variables for energies, number of basis functions and number of occupied orbitals
    E_HF, E_MP2, VNN, E_P, nbf, nocc = 0, 0, 0, 1, S.shape[0], sum(atoms) // 2

    # calculate nuclear-nuclear repulsion
    for i in range(len(atoms)):
        for j in range(i):
            VNN += atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :]) / A2BOHR

    # define Hamiltonian, matrix of coefficients, density matrix, exchange integral and orbital energies
    H, C, D, K, eps = T + V, np.zeros_like(S), np.zeros_like(S), J.transpose(0, 3, 2, 1), np.zeros([S.shape[0]])

    # while the difference in energy is greater than threshold
    while abs(E_HF - E_P) > thresh:
        # build the Fock matrix
        F = H + np.einsum("ijkl,ij", J - 0.5 * K, D, optimize=True)

        # solve the Fock equations
        eps, C = sp.linalg.eigh(F, S)

        # build the density from coefficients
        D = 2 * np.einsum("ij,kj", C[:, :nocc], C[:, :nocc])

        # calculate electron energy
        E_P, E_HF = E_HF, 0.5 * np.einsum("ij,ij", D, H + F, optimize=True)

    # transform the Hamiltonian and Coulomb integral to MO basis
    Hmo, Jmo = np.einsum('ip,ij,jq', C, H, C), np.einsum("ip,jq,ijkl,kr,ls", C, C, J, C, C, optimize=True)

    # calculate the MP2 correlation
    for i in range(nocc):
        for j in range(nocc):
            for a in range(nocc, len(H)):
                for b in range(nocc, len(H)):
                    E_MP2 += Jmo[i, a, j, b] * (2 * Jmo[i, a, j, b] - Jmo[i, b, j, a]) / (eps[i] + eps[j] - eps[a] - eps[b]);

    # # generate matrix of coefficients in MS basis
    # Cms = np.repeat(C, 2, axis=0); Cms = np.repeat(Cms, 2, axis=1)
    # # Cms = np.kron(np.eye(Cms.shape[0]), Cms)
    # spinind = np.arange(Cms.shape[0], dtype=int) % 2
    # print(np.repeat(np.repeat(np.kron(np.eye(C.shape[0]), np.ones(C.shape)), 2, axis=0), 2, axis=1))
    # Cms *= (spinind.reshape(-1, 1) == spinind)
    # print(Cms)
    #
    # # generate Hamiltonian in MS basis
    # Hms = np.repeat(Hmo, 2, axis=0); Hms = np.repeat(Hms, 2, axis=1)
    # spinind = np.arange(Hms.shape[0], dtype=int) % 2
    # Hms *= (spinind.reshape(-1, 1) == spinind)
    #
    # # generate coulomb integral in MS basis
    # # Tensor<4> Js = Eigen::Kron<4>(Matrix::Identity(2, 2), Eigen::Kron<4>(Matrix::Identity(2, 2), J).shuffle(Array<4>{3, 2, 1, 0}));
    # # test = 
    # # transform = np.kron(np.eye(2 * nocc), Jmo.reshape(2 * [Jmo.shape[0]**2]))
    # # print(Jmo.reshape(2 * [Jmo.shape[0]**2]))
    # print(np.kron(np.eye(2), np.kron(np.eye(2), J)).shape)
    # Jms = np.einsum("ip,jq,ijkl,kr,ls", Cms, Cms, np.kron(np.eye(2), np.kron(np.eye(2), J.transpose(3, 2, 1, 0).reshape(4, 4))).reshape([4, 4, 4, 4]), Cms, Cms, optimize=True)
    #
    # # determinant container
    # dets = list()
    #
    # # generate excitation determinants
    # for alpha in itertools.combinations(range(nbf), nocc):
    #     for beta in itertools.combinations(range(nbf), nocc):
    #         dets.append(Determinant(alpha, beta))
    #
    # # CI Hamiltonian
    # HCI = np.zeros([len(dets), len(dets)])
    #
    # # fill the hamiltonian
    # for i in range(len(dets)):
    #     for j in range(len(dets)):
    #         HCI[i, j] = dets[i].element(dets[j])
    #
    #
    # # print(Jms)

    # print the results
    print("HF ENERGY:", E_HF + VNN)
    print("MP2 ENERGY:", E_HF + E_MP2 + VNN)

#!/usr/bin/env python

import itertools as it, numpy as np

"""
This is an exercise script made for educational purposes. The script is supposed to calculate the Hartree-Fock energy and the MP2 correlation energy of a provided molecule.
The script expects that the molecule.xyz and all the integral files are present in the same directory as the script, but modification is trivial.
The script is divided into three parts. In the first part, the molecule and the integrals are loaded. This part should not be modified by a student.
In the second part, the Hartree-Fock calculation is performed. The student should fill the while loop with the correct calculations as well as the nuclear-nuclear repulsion energy.
In the last part, the MP2 calculation is performed. The student should transform the Coulomb integrals to the molecular orbital basis and calculate the MP2 correlation energy.
"""

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
    # OBTAIN THE MOLECULE AND ATOMIC INTEGRALS =========================================================================================================================================================

    # get the atomic numbers and coordinates of all atoms
    atoms = np.array([ATOM[line.split()[0]] for line in open("molecule.xyz").readlines()[2:]], dtype=int)
    coords = np.array([line.split()[1:] for line in open("molecule.xyz").readlines()[2:]], dtype=float)

    # convert to bohrs
    coords *= 1.8897261254578281

    # load the integrals from the files
    S, T, V = np.loadtxt("S.mat"), np.loadtxt("T.mat"), np.loadtxt("V.mat"); J = np.loadtxt("J.mat").reshape(4 * [S.shape[1]])

    # HARTREE-FOCK METHOD ==============================================================================================================================================================================

    """
    Here are some variables used throughout the calculation. Please use E_HF as a container for the HF energy. The E_HF_P is the previous energy and is used to check the convergence of the calculation.
    You can of course use different names for these variables, but some parts of the code are written with these names in mind. Threshold is the convergence threshold for the calculation.
    The nocc is the number of occupied orbitals for the system. The nbf is the number of basis functions. The variables E_HF and E_HF_P are initialized to zero and one, so that the SCF loop would start.
    """
    E_HF, E_HF_P, nocc, nbf, thresh = 0, 1, sum(atoms) // 2, S.shape[0], 1e-8

    """
    These lines define the Hamiltonian matrix as the sum of kinetic and potential matrix. The density matrix is initialized as a zero matrix and the coefficients are an empty array.
    The coefficient matrix is actually calculated in the while loop, but I already defined it here so it can be used outside of the loop (for example in the MP2 calculation).
    The exchange tensor is also correctly calculated here as a transpose of the Coulomb tensor. The eps vector is a vector of orbital energies and is defined here for the same reason as the coefficients.
    """
    K, eps = J.transpose(0, 3, 2, 1), np.array(nbf * [0])
    H, D, C = T + V, np.zeros_like(S), np.zeros_like(S)

    """
    This while loop is the SCF loop. Please fill it so it calculates the Fock matrix, solves the Fock equations, builds the density matrix from the coefficients and calculates the energy.
    You can use all the variables defined above and all the functions in numpy package. The recommended functions are np.einsum and np.linalg.eigh.
    Part of the calculation will probably be calculation of the inverse square root of a matrix. The numpy package does not conatin a function for this.
    You can find a library that can do that or you can do it manually. The manual calculation is, of course, preferred.
    """
    while abs(E_HF - E_HF_P) > thresh:
        break

    """
    In the followng block of code, please calculate the nuclear-nuclear repulsion energy. You should use only the atoms and coords variables.
    The code can be as short as two lines of. The result should be stored in the VNN variable.
    """
    VNN = 0

    # print the results
    print("RHF ENERGY: {:.8f}".format(E_HF + VNN))

    # MOLLER-PLESSET PERTRUBATION THEORY OF 2ND ORDER ==================================================================================================================================================

    """
    To perform the MP2 calculation, you will definitely need the Coulomb integrals in the molecular orbital basis. Please calculate them and store them in the Jmo variable.
    """
    Jmo = np.zeros_like(J)

    """
    And lastly, please calculate the MP2 correlation energy. The result should be stored in the E_MP2 variable.
    """
    E_MP2 = 0

    # print the results
    print("MP2 ENERGY: {:.8f}".format(E_HF + E_MP2 + VNN))

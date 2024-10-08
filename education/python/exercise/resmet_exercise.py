#!/usr/bin/env python

import itertools as it, numpy as np

"""
This educational script is designed to compute the Hartree-Fock energy, MP2 correlation energy and the FCI energy for a specified molecule. It requires that the molecule.xyz file and all integral files be located in the same directory as the script, though adjustments can easily be made. The script is structured into four main sections.

In the first section, the script loads the molecule and integral data. This section is intended to remain unaltered by students.

In the second section, the script performs the Hartree-Fock calculation. Students are expected to complete the while loop with the appropriate calculations. The calculation of the nuclear-nuclear repulsion energy is also requested.

The third section involves the transformation of the integrals to the molecular spinorbital basis. Here, students need to define the tiling matrix, spin masks, and transform the coefficient matrix to the molecular spinorbital basis. The Hamiltonian and Fock matrix in the molecular spinorbital basis should also be calculated.

The fourth section is dedicated to the MP2 calculation. Here, students need to transform the Coulomb integrals to the molecular orbital basis and compute the MP2 correlation energy.

The fifth section is dedicated to the CCSD calculation. Students should define the t1 and t2 amplitudes and complete the while loop with the appropriate calculations.

The sixth and final section involves the FCI (Full Configuration Interaction) calculation. Students should transform the required integrals into the molecular spinorbital basis, generate all possible electron determinants, and apply the Slater-Condon rules to construct the CI Hamiltonian.
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
    H, S = np.loadtxt("H_AO.mat", skiprows=1), np.loadtxt("S_AO.mat", skiprows=1); J = np.loadtxt("J_AO.mat", skiprows=1).reshape(4 * [S.shape[1]])

    # HARTREE-FOCK METHOD ==============================================================================================================================================================================

    """
    Here are defined some of the necessary variables. The variable "E_HF" stores the Hartree-Fock energy, while "E_HF_P" keeps track of the previous iteration's energy to monitor convergence. The "thresh" defines the convergence criteria for the calculation. The variables "nocc" and "nbf" represent the number of occupied orbitals and the number of basis functions, respectively. Initially, "E_HF" is set to zero and "E_HF_P" to one to trigger the start of the Self-Consistent Field (SCF) loop. Although you can rename these variables, it is important to note that certain sections of the code are tailored to these specific names.
    """
    E_HF, E_HF_P, nocc, nbf, thresh = 0, 1, sum(atoms) // 2, S.shape[0], 1e-8

    """
    These lines set up key components for our HF calculations. We initialize the density matrix as a zero matrix, and the coefficients start as an empty array. Although the coefficient matrix is computed within the while loop, it's defined outside to allow for its use in subsequent calculations, such as the MP energy computation. Similarly, the exchange tensor is accurately calculated here by transposing the Coulomb tensor. The "eps" vector, which contains the orbital energies, is also defined at this stage to facilitate access throughout the script. This setup ensures that all necessary variables are ready for iterative processing and further calculations beyond the SCF loop.
    """
    K, F, D, C, eps = J.transpose(0, 3, 2, 1), np.zeros_like(S), np.zeros_like(S), np.zeros_like(S), np.array(nbf * [0])

    """
    This while loop is the SCF loop. Please fill it so it calculates the Fock matrix, solves the Fock equations, builds the density matrix from the coefficients and calculates the energy. You can use all the variables defined above and all the functions in numpy package. The recommended functions are np.einsum and np.linalg.eigh. Part of the calculation will probably be calculation of the inverse square root of a matrix. The numpy package does not conatin a function for this. You can find a library that can do that or you can do it manually. The manual calculation is, of course, preferred.
    """
    while abs(E_HF - E_HF_P) > thresh:
        break

    """
    In the followng block of code, please calculate the nuclear-nuclear repulsion energy. You should use only the atoms and coords variables. The code can be as short as two lines. The result should be stored in the "VNN" variable.
    """
    VNN = 0

    # print the results
    print("RHF ENERGY: {:.8f}".format(E_HF + VNN))

    # INTEGRAL TRANSFORMS FOR POST-HARTREE-FOCK METHODS =================================================================================================================================================

    """
    To perform most of the post-HF calculations, we need to transform the Coulomb integrals to the molecular spinorbital basis, so if you don't plan to calculate any post-HF methods, you can end the eercise here. The restricted MP2 calculation could be done using the Coulomb integral in MO basis, but for the sake of subsequent calculations, we enforce here the integrals in the MS basis. The first thing you will need for the transform is the coefficient matrix in the molecular spinorbital basis. To perform this transform using the mathematical formulation presented in the materials, the first step is to form the tiling matrix "P" which will be used to duplicate columns of a general matrix. Please define it here.
    """
    P = np.zeros((nbf, 2 * nbf))

    """
    Now, please define the spin masks "M" and "N". These masks will be used to zero out spinorbitals, that should be empty.
    """
    M, N = np.zeros((nbf, 2 * nbf)), np.zeros((nbf, 2 * nbf))

    """
    With the tiling matrix and spin masks defined, please transform the coefficient matrix into the molecular spinorbital basis. The resulting matrix should be stored in the "Cms" variable.
    """
    Cms = np.zeros(2 * np.array(C.shape))

    """
    For some of the post-HF calculations, we will also need the Hamiltonian and Fock matrix in the molecular spinorbital basis. Please transform it and store it in the "Hms" and "Fms" variable. If you don't plan to calculate the CCSD method, you can skip the transformation of the Fock matrix, as it is not needed for the MP2 and CI calculations.
    """
    Hms, Fms = np.zeros(2 * np.array(H.shape)), np.zeros(2 * np.array(H.shape))

    """
    With the coefficient matrix in the molecular spinorbital basis available, we can proceed to transform the Coulomb integrals. It is important to note that the transformed integrals will contain twice as many elements along each axis compared to their counterparts in the atomic orbital (AO) basis. This increase is due to the representation of both spin states in the molecular spinorbital basis.
    """
    Jms = np.zeros(2 * np.array(J.shape))

    """
    The post-HF calculations also require the antisymmetrized two-electron integrals in the molecular spinorbital basis. These integrals are essential for the MP2 and CC calculations. Please define the "Jmsa" tensor as the antisymmetrized two-electron integrals in the molecular spinorbital basis.
    """
    Jmsa = np.zeros(2 * np.array(J.shape))

    """
    As mentioned in the materials, it is also practical to define the tensors of reciprocal orbital energy differences in the molecular spinorbital basis. These tensors are essential for the MP2 and CC calculations. Please define the "Emss", "Emsd" and "Emst" tensors as tensors of single, double and triple excitation energies, respectively. The configuration interaction will not need these tensors, so you can skip this step if you don't plan to program the CI method. The MP methods will require only the "Emsd" tensor, while the CC method will need both tensors.
    """
    Emss, Emsd = np.array([]), np.array([])

    # MOLLER-PLESSET PERTURBATION THEORY ===============================================================================================================================================================

    """
    Since we have everything we need for the MP calculations, we can now calculate the MP2 correlation energy. The result should be stored in the "E_MP2" variable.
    """
    E_MP2 = 0

    """
    Let's not stop here. We can calculate MP3 correlation energy as well. Please calculate it and store it in the "E_MP3" variable.
    """
    E_MP3 = 0

    # print the results
    print("MP2 ENERGY: {:.8f}".format(E_HF + E_MP2 +       + VNN))
    print("MP3 ENERGY: {:.8f}".format(E_HF + E_MP2 + E_MP3 + VNN))

    # COUPLED CLUSTER METHOD ===========================================================================================================================================================================

    """
    We also have everything we need for the CC calculations. In this exercise, we will calculate the CCSD energy. Since the calculation will be iterative, I define here the CCSD energy as zero, the "E_CCSD_P" variable will be used to monitor convergence.
    """
    E_CCSD, E_CCSD_P = 0, 1

    """
    The first step of the calculation is to define the "t1" and "t2" amplitudes. These arrays can be initialized as zero arrays with the appropriate dimensions. I will leave this task to you.
    """
    t1, t2 = np.array([]), np.array([])

    """
    Now for the more complicated part. The CCSD calculation is iterative, and the convergence criterion is set by the "thresh" variable. The while loop should be filled with the appropriate calculations. The calculation of the "t1" and "t2" amplitudes is the most challenging part of the CCSD calculation. After convergence, the "E_CCSD" variable should store the final CCSD energy.
    """
    while abs(E_CCSD - E_CCSD_P) > thresh:
        break

    # print the CCSD energy
    print("CCSD ENERGY: {:.8f}".format(E_HF + E_CCSD + VNN))

    # CONFIGURATION INTERACTION ========================================================================================================================================================================

    """
    Since we already calculated the necessary integrals in the MS basis, we can proceed. The next step involves generating determinants. We will store these in a simple list, with each determinant represented by an array of numbers, where each number corresponds to an occupied spinorbital. Since we are programming for Full Configuration Interaction (FCI), we aim to generate all possible determinants. However, should we decide to implement methods like CIS, CID, or CISD, we could easily limit the number of excitations. It is important to remember that for all CI methods, the rest of the code remains unchanged. The only difference lies in the determinants used. Don't overcomplicate this. Generating all possible determinants can be efficiently achieved using a simple list comprehension. I recommend employing the combinations function from the itertools package to facilitate this task.
    """
    dets = list()

    """
    Now, for your convenince, I define here the CI Hamiltonian.
    """
    Hci = np.zeros([len(dets), len(dets)])

    """
    Before we begin constructing the Hamiltonian, I recommend defining the Slater-Condon rules. Let's consider that the input for these functions will be an array of spinorbitals, segmented into unique and common ones. A practical approach might be to arrange this 1D array with all unique spinorbitals at the front, followed by the common spinorbitals. This arrangement allows you to easily determine the number of unique spinorbitals based on the rule being applied, meaning you will always know how many entries at the beginning of the array are unique spinorbitals. While you can develop your own method for managing this array, I will proceed under the assumption that the Slater-Condon rules we use will take a single array of spinorbitals and return an unsigned matrix element. The sign of this element will be corrected later in the script. For simplicity and flexibility, I'll define these rules using lambda functions, but you're welcome to expand them into full functions if you prefer.
    """
    slater0 = lambda so: 0
    slater1 = lambda so: 0
    slater2 = lambda so: 0

    """
    We can now proceed to filling the CI Hamiltonian. The loop is simple.
    """
    for i in range(Hci.shape[0]):
        for j in range(Hci.shape[1]):

            """
            The challenging part of this process is aligning the determinants. In this step, I transfer the contents of the j-th determinant into the "aligned" determinant. It's important not to alter the j-th determinant directly within its original place, as doing so could disrupt the computation of other matrix elements. Instead, we carry out the next steps on the determinant now contained in the "aligned" variable. Additionally, the element sign is defined at this stage. You probably want to leave this unchanged.
            """
            aligned, sign = dets[j].copy(), 1

            """
            Now it's your turn. Please adjust the "aligned" determinant to match the i-th determinant as closely as possible. By "align", I mean you should execute a series of spinorbital swaps to minimize the differences between the "aligned" and the i-th determinant. It's also important to monitor the number of swaps you make, as each swap affects the sign of the determinant, hence the reason for the "sign" variable defined earlier. This task is not straightforward, so don't hesitate to reach out to the authors if you need guidance.
            """
            aligned = aligned

            """
            After aligning, we end up with two matched determinants: "aligned" and "dets[i]". At this point, we can apply the Slater-Condon rules. I suggested earlier that the input for these rules should be an array combining both unique and common spinorbitals. You can prepare this array now. However, if you've designed your Slater-Condon rules to directly accept the determinants instead, you can skip this preparatory step.
            """
            so = list()

            """
            Now, you'll need to assign the matrix element. Start by determining the number of differences between the two determinants. Based on this number, apply the corresponding Slater-Condon rule. Don't forget to multiply the result by the sign to account for any changes due to swaps made during the alignment of the determinants.
            """
            H[i, j] = 0

    """
    You can finally solve the eigenvalue problem. Please, assign the correlation energy to the "E_FCI" variable.
    """
    E_FCI = 0

    # print the results
    print("FCI ENERGY: {:.8f}".format(E_HF + E_FCI + VNN))

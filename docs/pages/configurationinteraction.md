---
title: Configuration Interaction
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Configuration Interaction

Configuration Interaction is a post-Hartree--Fock, utilizing a linear variational approach to address the nonrelativistic Schrödinger equation under the Born–Oppenheimer approximation for multi-electron quantum systems. Configuration Interaction mathematically represents the wave function as a linear combination of Slater determinants. The term "configuration" refers to different ways electrons can occupy orbitals, while "interaction" denotes the mixing of these electronic configurations or states. Configuration Interaction computations, however, are resource-intensive, requiring significant CPU time and memory, limiting their application to smaller molecular systems. While Full Configuration Interaction considers all possible electronic configurations, making it computationally prohibitive for larger systems, truncated versions like Configuration Interaction Singles and Doubles or Configuration Interaction Singles, Doubles and Triples are more feasible and commonly employed in quantum chemistry studies.

## Theoretical Background of General Configuration Interaction

The idea is quite simple, using the convention, that the indices $i$, $j$, $k$, and $l$ run over occupied spinorbitals and the indices $a$, $b$, $c$, and $d$ run over virtual spinorbitals. The Configuration Interaction wavefunction is written as

\begin{equation}
\ket{\Psi}=c\_0\ket{\Psi\_0}+\left(\frac{1}{1!}\right)^2c\_i^a\ket{\Psi\_i^a}+\left(\frac{1}{2!}\right)^2c\_{ij}^{ab}\ket{\Psi\_{ij}^{ab}}+\left(\frac{1}{3!}\right)^2c\_{ijk}^{abc}\ket{\Psi\_{ijk}^{abc}}+\dots
\end{equation}

and we would like to know the coefficients $c$ that minimize the energy. To do that, we simply construct the hamiltonian matrix in the basis of excited determinants and diagonalize it. The Configuration Interaction Hamiltonian matrix $\mathbf{H}^{CI}$ is constructed as

\begin{equation}\label{eq:ci-hamiltonian}
\mathbf{H}^{CI}=
\begin{bmatrix}
\braket{\Psi\_0|\hat{H}|\Psi\_0} & \braket{\Psi\_0|\hat{H}|\Psi\_i^a} & \braket{\Psi\_0|\hat{H}|\Psi\_{ij}^{ab}} & \dots \\\\\
\braket{\Psi\_i^a|\hat{H}|\Psi\_0} & \braket{\Psi\_i^a|\hat{H}|\Psi\_i^a} & \braket{\Psi\_i^a|\hat{H}|\Psi\_{ij}^{ab}} & \dots \\\\\
\braket{\Psi\_{ia}^{jb}|\hat{H}|\Psi\_0} & \braket{\Psi\_{ia}^{jb}|\hat{H}|\Psi\_1} & \braket{\Psi\_{ia}^{jb}|\hat{H}|\Psi\_{ia}^{jb}} & \dots \\\\\
\vdots & \vdots & \vdots & \ddots
\end{bmatrix}
\end{equation}

and solving the eigenvalue problem

\begin{equation}\label{eq:ci-eigenvalue-problem}
\mathbf{H}^{CI}\mathbf{C}^{CI}=\mathbf{C}^{CI}\mathbf{\varepsilon}^{CI}
\end{equation}

where $\mathbf{C}^{CI}$ is a matrix of coefficients and $\mathbf{\varepsilon}^{CI}$ is a diagonal matrix of eigenvalues. The lowest eigenvalue corresponds to the ground state energy, while the eigenvector corresponding to the lowest eigenvalue gives the coefficients that minimize the energy. The matrix elements of the Configuration Interaction Hamiltonian are calculated using the Slater-Condon rules in the form

\begin{equation}\label{eq:slater-condon-rules}
\mathbf{H}_{ij}^{CI}=
\begin{cases} 
\displaystyle \sum_kH\_{kk}^{core,MS}+\frac{1}{2}\sum_k\sum_l\braket{kl||kl}&D_i=D_j \\\\\
\displaystyle H\_{pr}^{core,MS}+\sum_k\braket{pk||lk}&D_i=\left\lbrace\dotsi p\dotsi\right\rbrace\land D_j=\left\lbrace\dotsi r\dotsi\right\rbrace \\\\\
\displaystyle \vphantom{\sum_k}\braket{pq||rs}&D_i=\left\lbrace\dotsi p\dotsi q\dotsi\right\rbrace\land D_j=\left\lbrace\dotsi r\dotsi s\dotsi\right\rbrace \\\\\
\displaystyle \vphantom{\sum_k}0&\text{otherwise},
\end{cases}
\end{equation}

where $D\_i$ and $D\_j$ are Slater determinants, $\mathbf{H}^{core,MS}$ is the core Hamiltonian in the basis of molecular spinorbitals, and $\braket{pk\|\|lk}$ are the antisymmetrized Coulomb repulsion integrals in the basis of molecular spinorbitals and physicists' notation. The sums extend over all spinorbitals common between the two determinants. All the integrals in the MS basis are already explained [here](hartreefockmethod.html#integral-transforms-to-the-basis-of-molecular-spinorbitals). Keep in mind, that to apply the Slater-Condon rules, the determinants must be aligned, and the sign of the matrix elements must be adjusted accordingly, based on the number of permutations needed to align the determinants.

The problem with Configuration Interaction is that it is not size-extensive, meaning that the energy does not scale linearly with the number of electrons. This is because the Configuration Interaction wavefunction is not size-consistent, and the energy of a system is not the sum of the energies of its parts. This is a significant drawback of Configuration Interaction, as it limits its application to small systems.

## Full Configuration Interaction Implementation

Let's consider the Full Configuration Interaction method, which considers all possible electronic configurations within a given basis set. The Full Configuration Interaction method provides the most accurate description of the electronic structure, but its computational cost grows exponentially with the number of electrons and basis functions, making it infeasible for large systems.

The only thing that is needed, besides the general Configuration Interaction equations, are the determinants. For simplicity, we will include singlet and triplet states. The number of these determinants $N\_D$ can be for Full Configuration Interaction calculated using the binomial coefficients

\begin{equation}
N\_D=\binom{n}{k}
\end{equation}

assuming $k$ is the total number of electrons, and $n$ is the total number of spinorbitals. Each determinant is formed by permuting the electrons between spinorbitals. For practical representation, it's useful to describe determinants as arrays of numbers, where each number corresponds to the index of an occupied orbitals. For example, the ground state determinant for a system with 6 electrons and 10 spinorbitals can be represented as $\left\lbrace 0,1,2,3,4,5\right\rbrace$, whereas the determinant $\left\lbrace 0,1,2,3,4,6\right\rbrace$ represents an excited state with one electron excited from orbital 5 to orbital 6. Using the determinants, the Configuration Interaction Hamiltonian matrix \eqref{eq:ci-hamiltonian} can be constructed, and the eigenvalue problem \eqref{eq:ci-eigenvalue-problem} can be solved to obtain the ground and excited state energies.

## Full Configuration Interaction Code Exercise

This section provides a simple example of a Full Configuration Interaction implementation. The code snippets are, as the previous coding examples, written in Python and use the NumPy library for numerical operations. The example demonstrates the generation of determinants, the construction of the Configuration Interaction Hamiltonian matrix, and the solution of the eigenvalue problem to obtain the ground state energy. The code is designed for educational purposes and may require modifications for practical applications, such as the inclusion of truncation schemes for larger systems. You are expected to have the Hartree--Fock calculations [here](hartreefockmethod.html#hartreefock-method-and-integral-transform-coding-exercise) completed before proceeding with the Configuration Interaction implementation, as the Configuration Interaction method builds on the Hartree--Fock method.

In this exercise, you are expected to calculate the Full Configuration Interaction energy for a simple system. The exercise is provided in the Listing <!--\ref{code:ci_exercise}--> below.

<!--{id=code:ci_exercise caption="Configuration Interaction exercise code."}-->
```python
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
```

Solution to this exercise can be found [here](codesolutions.html#full-configuration-interaction).

{:.cite}

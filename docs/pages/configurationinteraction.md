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

The idea is quite simple, using the convention, that the indices $i$, $j$, $k$, and $l$ run over occupied spinorbitals and the indices $a$, $b$, $c$, and $d$ run over virtual spinorbitals. The CI wavefunction is written as

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

where $\mathbf{C}^{CI}$ is a matrix of coefficients and $\mathbf{\varepsilon}^{CI}$ is a diagonal matrix of eigenvalues. The lowest eigenvalue corresponds to the ground state energy, while the eigenvector corresponding to the lowest eigenvalue gives the coefficients that minimize the energy. The matrix elements of the CI Hamiltonian are calculated using the Slater-Condon rules in the form

\begin{equation}\label{eq:slater-condon-rules}
\mathbf{H}_{ij}^{CI}=
\begin{cases} 
\displaystyle \sum_kH\_{kk}^{core,MS}+\frac{1}{2}\sum_k\sum_l\braket{kl||kl}&D_i=D_j \\\\\
\displaystyle H\_{pr}^{core,MS}+\sum_k\braket{pk||lk}&D_i=\left\lbrace\dotsi p\dotsi\right\rbrace\land D_j=\left\lbrace\dotsi r\dotsi\right\rbrace \\\\\
\displaystyle \vphantom{\sum_k}\braket{pq||rs}&D_i=\left\lbrace\dotsi p\dotsi q\dotsi\right\rbrace\land D_j=\left\lbrace\dotsi r\dotsi s\dotsi\right\rbrace \\\\\
\displaystyle \vphantom{\sum_k}0&\text{otherwise},
\end{cases}
\end{equation}

where $D\_i$ and $D\_j$ are Slater determinants, $\mathbf{H}^{core,MS}$ is the core Hamiltonian in the basis of molecular spinorbitals, and $\braket{pk\|\|lk}$ are the antisymmetrized Coulomb repulsion integrals in the basis of molecular spinorbitals and physicists' notation. The sums extend over all spinorbitals common between the two determinants. All the integrals in the MS basis are already explained [here](hartreefockmethod.html#the-integral-transforms). Keep in mind, that to apply the Slater-Condon rules, the determinants must be aligned, and the sign of the matrix elements must be adjusted accordingly, based on the number of permutations needed to align the determinants.

The problem with Configuration Interaction is that it is not size-extensive, meaning that the energy does not scale linearly with the number of electrons. This is because the Configuration Interaction wavefunction is not size-consistent, and the energy of a system is not the sum of the energies of its parts. This is a significant drawback of Configuration Interaction, as it limits its application to small systems.

## Full Configuration Interaction Implementation

Let's consider the Full Configuration Interaction method, which considers all possible electronic configurations within a given basis set. The Full Configuration Interaction method provides the most accurate description of the electronic structure, but its computational cost grows exponentially with the number of electrons and basis functions, making it infeasible for large systems.

The only thing that is needed, besides the general Configuration Interaction equations, are the determinants. For simplicity, we will include singlet and triplet states. The number of these determinants $N\_D$ can be for Full Configuration Interaction calculated using the binomial coefficients

\begin{equation}
N\_D=\binom{n}{k}
\end{equation}

assuming $k$ is the total number of electrons, and $n$ is the total number of spinorbitals. Each determinant is formed by permuting the electrons between spinorbitals. For practical representation, it's useful to describe determinants as arrays of numbers, where each number corresponds to the index of an occupied orbitals. For example, the ground state determinant for a system with 6 electrons and 10 spinorbitals can be represented as $\left\lbrace 0,1,2,3,4,5\right\rbrace$, whereas the determinant $\left\lbrace 0,1,2,3,4,6\right\rbrace$ represents an excited state with one electron excited from orbital 5 to orbital 6. Using the determinants, the Configuration Interaction Hamiltonian matrix \ref{eq:ci-hamiltonian} can be constructed, and the eigenvalue problem \ref{eq:ci-eigenvalue-problem} can be solved to obtain the ground and excited state energies.
{:.cite}

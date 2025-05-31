---
title: Configuration Interaction
parent: Electronic Structure Methods
layout: default
nav_order: 3
---
{% include mathjax.html %}

# Configuration Interaction

Configuration Interaction is a post-Hartree--Fock, utilizing a linear variational approach to address the nonrelativistic Schrödinger equation under the Born–Oppenheimer approximation for multi-electron quantum systems. Configuration Interaction mathematically represents the wave function as a linear combination of Slater determinants. The term "configuration" refers to different ways electrons can occupy orbitals, while "interaction" denotes the mixing of these electronic configurations or states. Configuration Interaction computations, however, are resource-intensive, requiring significant CPU time and memory, limiting their application to smaller molecular systems. While Full Configuration Interaction considers all possible electronic configurations, making it computationally prohibitive for larger systems, truncated versions like Configuration Interaction Singles and Doubles or Configuration Interaction Singles, Doubles and Triples are more feasible and commonly employed in quantum chemistry studies.

## Theoretical Background of General Configuration Interaction

In Configuration Interaction theory, we expand the wavefunction $$\ket{\Psi}$$ in terms of the Hartree--Fock reference determinant and its excited configurations as

$$
\begin{equation}
\ket{\Psi}=c_0\ket{\Psi_0}+\left(\frac{1}{1!}\right)^2c_i^a\ket{\Psi_i^a}+\left(\frac{1}{2!}\right)^2c_{ij}^{ab}\ket{\Psi_{ij}^{ab}}+\left(\frac{1}{3!}\right)^2c_{ijk}^{abc}\ket{\Psi_{ijk}^{abc}}+\dots
\end{equation}
$$

where we seek the coefficients $$\mathbf{c}$$ that minimize the energy. To determine these coefficients, we construct and diagonalize the Hamiltonian matrix in the basis of these excited determinants. The Configuration Interaction Hamiltonian matrix $$\mathbf{H}^{\mathrm{CI}}$$ is represented as

$$
\begin{equation}\label{eq:ci-hamiltonian}
\mathbf{H}^{\mathrm{CI}}=
\begin{bmatrix}
\braket{\Psi_0\vert\hat{H}\vert\Psi_0} & \braket{\Psi_0\vert\hat{H}\vert\Psi_i^a} & \braket{\Psi_0\vert\hat{H}\vert\Psi_{ij}^{ab}} & \dots \\
\braket{\Psi_i^a\vert\hat{H}\vert\Psi_0} & \braket{\Psi_i^a\vert\hat{H}\vert\Psi_i^a} & \braket{\Psi_i^a\vert\hat{H}\vert\Psi_{ij}^{ab}} & \dots \\
\braket{\Psi_{ia}^{jb}\vert\hat{H}\vert\Psi_0} & \braket{\Psi_{ia}^{jb}\vert\hat{H}\vert\Psi_1} & \braket{\Psi_{ia}^{jb}\vert\hat{H}\vert\Psi_{ia}^{jb}} & \dots \\
\vdots & \vdots & \vdots & \ddots
\end{bmatrix}
\end{equation}
$$

After the Hamiltonian matrix is constructed, we solve the eigenvalue problem

$$
\begin{equation}\label{eq:ci-eigenvalue-problem}
\mathbf{H}^{\mathrm{CI}}\mathbf{C}^{\mathrm{CI}}=\mathbf{C}^{\mathrm{CI}}\mathbf{\varepsilon}^{\mathrm{CI}}
\end{equation}
$$

where $$\mathbf{C}^{\mathrm{CI}}$$ is a matrix of coefficients and $$\mathbf{\varepsilon}^{\mathrm{CI}}$$ is a diagonal matrix of eigenvalues. The lowest eigenvalue gives the ground-state energy, and the corresponding eigenvector provides the coefficients that minimize the energy. The elements of the Configuration Interaction Hamiltonian matrix are computed using the Slater--Condon rules, summarized in one function as

$$
\begin{equation}\label{eq:slater-condon-rules}
\mathbf{H}_{ij}^{\mathrm{CI}}=
\begin{cases} 
\displaystyle \sum_kH_{kk}^{\mathrm{core},\mathrm{MS}}+\frac{1}{2}\sum_k\sum_l\braket{kl\vert\vert kl}&D_i=D_j \\
\displaystyle H_{pr}^{\mathrm{core},\mathrm{MS}}+\sum_k\braket{pk\vert\vert rk}&D_i=\left\lbrace\dotsi p\dotsi\right\rbrace\land D_j=\left\lbrace\dotsi r\dotsi\right\rbrace \\
\displaystyle \vphantom{\sum_k}\braket{pq\vert\vert rs}&D_i=\left\lbrace\dotsi p\dotsi q\dotsi\right\rbrace\land D_j=\left\lbrace\dotsi r\dotsi s\dotsi\right\rbrace \\
\displaystyle \vphantom{\sum_k}0&\text{otherwise}
\end{cases}
\end{equation}
$$

where $$D_i$$ and $$D_j$$ are Slater determinants, $$\mathbf{H}^{\mathrm{core},\mathrm{MS}}$$ is the core Hamiltonian in the Molecular Spinorbital basis, and $$\braket{pk\vert\vert lk}$$ are the antisymmetrized two-electron integrals in Molecular Spinorbital basis and physicists' notation. The sums extend over all spinorbitals common between the two determinants. These integrals were previously transformed [here](hartreefock.html#integral-transforms-to-the-basis-of-molecular-spinorbitals)<!--in Section \ref{sec:integral_transform}-->. Keep in mind, that to apply the Slater-Condon rules, the determinants must be aligned, and the sign of the matrix elements must be adjusted accordingly, based on the number of permutations needed to align the determinants.

An important caveat in Configuration Interaction theory is its lack of size-extensivity, which implies that the energy does not scale linearly with the number of electrons. This drawback stems from the fact that the Configuration Interaction wavefunction is not size-consistent, meaning the energy of a combined system is not simply the sum of the energies of its isolated parts. This limitation restricts the application of Configuration Interaction mainly to small molecular systems.

## Full Configuration Interaction Implementation

In Full Configuration Interaction, we aim to account for all possible electronic configurations within a chosen basis set, offering the most accurate wavefunction representation for the given basis. Although this method yields highly precise electronic structure information, it is computationally intensive. Its cost scales exponentially with both the number of electrons and basis functions, limiting its feasibility to smaller systems.

The Full Configuration Interaction process involves constructing all possible Slater determinants for a system. For simplicity, we’ll assume that we want to include both singlet and triplet states in our determinant space. The total number of these determinants $$N_D$$ can be calculated using binomial coefficients

$$
\begin{equation}
N_D=\binom{n}{k}
\end{equation}
$$

where $$k$$ is the total number of electrons, and $$n$$ is the total number of spinorbitals. For practical representation, it's useful to describe determinants as arrays of numbers, where each number corresponds to the index of an occupied orbitals. For example, the ground state determinant for a system with 6 electrons can be represented as $$\left\lbrace 0,1,2,3,4,5\right\rbrace$$, whereas the determinant $$\left\lbrace 0,1,2,3,4,6\right\rbrace$$ represents an excited state with one electron excited from orbital 5 to orbital 6. Using the determinants, the Configuration Interaction Hamiltonian matrix \eqref{eq:ci-hamiltonian} can be constructed, and the eigenvalue problem \eqref{eq:ci-eigenvalue-problem} can be solved to obtain the ground and excited state energies.

{:.cite}

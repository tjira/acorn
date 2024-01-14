---
title: Hartree–Fock Method
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Hartree-Fock Method

The Hartree-Fock method stands as a cornerstone in quantum chemistry, offering a systematic approach to solve the electronic structure problem in molecules. This computational technique strives to determine the optimal wave function for a given molecular system, providing insights into the distribution of electrons and their energies. The method's foundation lies in the assumption of a single determinantal wave function, which embodies the antisymmetry principle inherent in quantum mechanics.

## Restricted Hartree–Fock Method

The restricted version of the Hartree-Fock method (RHF) introduces a constraint on electron spin, limiting the system to a fixed number of alpha and beta electrons. This approach ensures that the overall electron configuration adheres to the Pauli exclusion principle, preventing two electrons with the same spin from occupying the same orbital. While this simplification reduces computational complexity, it may not capture certain electronic correlations and is particularly suitable for closed-shell systems.

### The Method

Let's begin by defining the core Hamiltonian, also known as the one-electron Hamiltonian. The core Hamiltonian represents a part of the full Hamiltonian that excludes electron-electron repulsion. In index notation, it is expressed as 

\begin{equation}
H_{\mu\nu}^{core}=T_{\mu\nu}+V_{\mu\nu},
\end{equation}

where $\mu$ and $\nu$ are indices of atomic orbitals, $T_{\mu\nu}$ is a kinetic energy matrix element and $V_{\mu\nu}$ is a potential energy matrix element. These matrix elements, given in bra-ket notation as 

\begin{align}
T_{\mu\nu}&=\braket{\phi_{\mu}|\hat{T}|\phi_{\nu}} \\\\\
V_{\mu\nu}&=\braket{\phi_{\mu}|\hat{V}|\phi_{\nu}}
\end{align}

are calculated using a basis set and handled by libint library. Additionally, the overlap matrix elements

\begin{equation}
S_{\mu\nu}=\braket{\phi_{\mu}|\phi_{\nu}}
\end{equation}

and the two-electron Coulomb repulsion integral

\begin{equation}
J_{\mu\nu\kappa\lambda}=\braket{\phi_{\mu}\phi_{\mu}|\hat{J}|\phi_{\kappa}\phi_{\lambda}}
\end{equation}

play crucial roles. The Hartree-Fock method revolves around solving the Roothaan equations

\begin{equation}
\mathbf{FC}=\mathbf{SC\varepsilon},
\end{equation}

where $\mathbf{F}$ is the Fock matrix, $\mathbf{C}$ is a matrix of orbital coefficients and $\mathbf{\varepsilon}$ represents orbital energies. The problem is that the Fock matrix in the form

\begin{equation}
F_{\mu\nu}=H_{\mu\nu}^{core}+D_{\kappa\lambda}(J_{\mu\nu\kappa\lambda}-\frac{1}{2}J_{\mu\lambda\kappa\nu})
\end{equation}

However, the challenge lies in the iterative nature of the solution, as the Fock matrix depends on the unknown density matrix $D$. This iterative process is carried out through a self-consistent field (SCF) method. In each iteration, a guess for the density matrix is made, and the Roothaan equations are solved. The updated density matrix

\begin{equation}
D_{\mu\nu}=2C_{\mu i}C_{\nu i}
\end{equation}

and the total energy of the system

\begin{equation}
E=\frac{1}{2}D_{\mu\nu}(H_{\mu\nu}^{core}+F_{\mu\nu})+E_{nuc},
\end{equation}

are then calculated using the core Hamiltonian and the Fock matrix. The nuclear repulsion energy

\begin{equation}
E_{nuc}=\sum_{A}\sum_{B<A}\frac{Z_{A}Z_{B}}{R_{AB}}.
\end{equation}

is precomputed before the SCF iterations. The SCF procedure continues iteratively until a convergence criterion, often based on the difference in energy between consecutive iterations or the Frobenius norm of the density matrix, is met.

### The Gradient

If we perform the calculation as described above and get the density matrix $\mathbf{D}$ we can evaluate the nuclear energy gradient as

\begin{equation}
\frac{\partial E}{\partial X_{A,1}}=D_{\mu\nu}\frac{\partial H_{\mu\nu}^{core}}{\partial X_{A,1}}+2D_{\mu\nu}D_{\kappa\lambda}\frac{\partial J_{\mu\nu\kappa\lambda}}{\partial X_{A,1}}-2W_{\mu\nu}\frac{\partial S_{\mu\nu}}{\partial X_{A,1}},
\end{equation}

where $\mathbf{W}$ is energy weighed density matrix defined as

\begin{equation}
W_{\mu\nu}=2C_{\mu i}C_{\nu i}\varepsilon_i.
\end{equation}

## Unrestricted Hartree–Fock Method

In contrast to the RHF, the unrestricted version of the Hartree-Fock method removes the spin restriction, allowing alpha and beta electrons to occupy the same set of molecular orbitals. This extension enables a more accurate representation of open-shell systems, where the number of electrons with opposite spins may differ. While the unrestricted approach introduces computational challenges, it proves invaluable in describing molecules with unpaired electrons, providing a comprehensive understanding of their electronic structure and properties.

### The Method

As the one and two-electron integrals remain unaffected by the density matrix, they retain their form from the restricted version. We introduce the $\alpha$ and $\beta$ Fock matrices defined as

\begin{align}
F_{\mu\nu}^{\alpha}=H_{\mu\nu}^{core}+\frac{1}{2}D_{\kappa\lambda}^{\alpha}J_{\mu\nu\kappa\lambda}+\frac{1}{2}D_{\kappa\lambda}^{\beta}J_{\mu\nu\kappa\lambda}-\frac{1}{2}D_{\kappa\lambda}^{\alpha}J_{\mu\lambda\kappa\nu} \\\\\
F_{\mu\nu}^{\beta}=H_{\mu\nu}^{core}+\frac{1}{2}D_{\kappa\lambda}^{\alpha}J_{\mu\nu\kappa\lambda}+\frac{1}{2}D_{\kappa\lambda}^{\beta}J_{\mu\nu\kappa\lambda}-\frac{1}{2}D_{\kappa\lambda}^{\beta}J_{\mu\lambda\kappa\nu}
\end{align}

are formed similarly to the restricted version, with a small distinction. It is crucial to handle matrix slicing with care, ensuring that the number of columns in these matrices matches the number of corresponding electrons. Using these density matrices, we solve two distinct Fock equations, acquiring expansion coefficients and orbital energies before updating the density matrices. The total electron energy is computed as

\begin{equation}
E=\frac{1}{4}D_{\mu\nu}^{\alpha}(H_{\mu\nu}^{core}+F_{\mu\nu}^{\alpha})+\frac{1}{4}D_{\mu\nu}^{\beta}(H_{\mu\nu}^{core}+F_{\mu\nu}^{\beta})+E_{nuc}.
\end{equation}

This formulation maintains the same nuclear energy term as in the restricted version.

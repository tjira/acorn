---
title: Hartree–Fock Method
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Hartree--Fock Method

The Hartree--Fock method stands as a cornerstone in quantum chemistry, offering a systematic approach to solve the electronic structure problem in molecules. This computational technique strives to determine the optimal wave function for a given molecular system, providing insights into the distribution of electrons and their energies.

## Theoretical Background

Ultimately, we are interested in solving the Schrödinger equation in the form

\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi},
\end{equation}

where $\hat{H}$ is the molecular Hamiltonian operator, $\Psi$ is the molecular wave function, and $E$ is the total energy of the system. The Hartree--Fock method aims to approximate the wave function $\Psi$ by a single Slater determinant, which we will write in the form

\begin{equation}
\ket{\Psi}=\ket{\chi\_1\chi\_2\cdots\chi\_N},
\end{equation}

where $\chi\_i$ represents a Molecular Spinorbital and N is the total number of electrons. The Hartree--Fock method seeks to optimize these molecular orbitals to minimize the total energy of the system, providing a reliable estimate of the electronic structure. However, in the Restricted Hartree--Fock method, we impose a constraint on the electron spin, and instead of spin orbitals, we work with spatial orbitals, which allows us to write the slater determinant in terms of spatial orbitals in the form

\begin{equation}
\ket{\Psi}=\ket{\Phi\_1\Phi\_2\cdots\Phi\_{N/2}},
\end{equation}

where $\Phi\_i$ represents a molecular spatial orbital. We can see, that for the Restricted Hartree--Fock method, we need to have an even number of electrons in the system. In practical calculations, it is convenient to expand the molecular orbitals (spin or spatial) in terms of basis functions $\phi$ (usually sums of Gaussian functions) and work with the expansion coefficients. If we assume that the wavefunction is a single Slater determinant, our molecular orbitals are expanded in terms of basis functions and optimize the energy of such determinant, we arrive at the Roothaan equations in the form

\begin{equation}\label{eq:roothaan}
\mathbf{FC}=\mathbf{SC\varepsilon},
\end{equation}

where $\mathbf{F}$ is the Fock matrix (defined later), $\mathbf{C}$ is a matrix of orbital coefficients, $\mathbf{S}$ is the overlap matrix (also defined later), and $\mathbf{\varepsilon}$ represents orbital energies.

## Implementation of the Restricted Hartree–Fock Method

Let's begin by defining the core Hamiltonian, also known as the one-electron Hamiltonian. The core Hamiltonian represents a part of the full Hamiltonian that excludes electron-electron repulsion. In index notation, it is expressed as 

\begin{equation}\label{eq:hamiltonian}
H\_{\mu\nu}^{core}=T\_{\mu\nu}+V\_{\mu\nu}
\end{equation}

where $\mu$ and $\nu$ are indices of basis functions, $T\_{\mu\nu}$ is a kinetic energy matrix element and $V\_{\mu\nu}$ is a potential energy matrix element. These matrix elements are given as

\begin{align}
T\_{\mu\nu}&=\braket{\phi\_{\mu}|\hat{T}|\phi\_{\nu}} \\\\\
V\_{\mu\nu}&=\braket{\phi\_{\mu}|\hat{V}|\phi\_{\nu}}
\end{align}

are usually calculated using analytical expressions. Additionally, using analytical expressions, we can calculate the overlap integrals

\begin{equation}\label{eq:overlap}
S\_{\mu\nu}=\braket{\phi\_{\mu}|\phi\_{\nu}}
\end{equation}

and the two-electron Coulomb repulsion integrals

\begin{equation}\label{eq:coulomb}
J\_{\mu\nu\kappa\lambda}=\braket{\phi\_{\mu}\phi\_{\mu}|\hat{J}|\phi\_{\kappa}\phi\_{\lambda}},
\end{equation}

which play crucial roles in the Hartree--Fock calculation.<!--\cite{10.1016/S0065-3276!08!60019-2}--> The Hartree--Fock method revolves around solving the Roothaan equations \ref{eq:roothaan} iteratively. The Fock matrix is defined as

\begin{equation}\label{eq:fock}
F\_{\mu\nu}=H\_{\mu\nu}^{core}+D\_{\kappa\lambda}(J\_{\mu\nu\kappa\lambda}-\frac{1}{2}J\_{\mu\lambda\kappa\nu})
\end{equation}

depends on the unknown density matrix $\mathbf{D}$. This iterative process is carried out through a Self-Consistent Field method. In each iteration, a guess for the density matrix is made, and the Roothaan equations are solved. The density matrix is then updated as

\begin{equation}
D\_{\mu\nu}=2C\_{\mu i}C\_{\nu i}
\end{equation}

and the total energy of the system

\begin{equation}
E=\frac{1}{2}D\_{\mu\nu}(H\_{\mu\nu}^{core}+F\_{\mu\nu})+E\_{nuc}
\end{equation}

is then calculated using the core Hamiltonian and the Fock matrix. The important thing to note is that the Fock matrix depends on the density matrix, which is updated in each iteration. The initial guess for the density matrix is often set to zero. After the density matrix and the total energy are converged, the Self-Consistent Field procedure is terminated, and the optimized molecular orbitals are obtained. To get the final energy of the system, we need to add the nuclear repulsion energy of the form

\begin{equation}
E\_{nuc}=\sum\_{A}\sum\_{B<A}\frac{Z\_{A}Z\_{B}}{R\_{AB}},
\end{equation}

where $Z\_A$ is the nuclear charge of atom A, and $R\_{AB}$ is the distance between atoms A and B.

### Gradient of the Restricted Hartree--Fock Method

If we perform the calculation as described above and get the density matrix $\mathbf{D}$ we can evaluate the nuclear energy gradient as

\begin{equation}
\frac{\partial E}{\partial X\_{A,i}}=D\_{\mu\nu}\frac{\partial H\_{\mu\nu}^{core}}{\partial X\_{A,i}}+2D\_{\mu\nu}D\_{\kappa\lambda}\frac{\partial J\_{\mu\nu\kappa\lambda}}{\partial X\_{A,i}}-2W\_{\mu\nu}\frac{\partial S\_{\mu\nu}}{\partial X\_{A,i}}
\end{equation}

where $i$ is the index of the coordinate and where $\mathbf{W}$ is energy weighed density matrix defined as

\begin{equation}
W_{\mu\nu}=2C_{\mu i}C_{\nu i}\varepsilon_i
\end{equation}

## Integral Transforms to the Basis of Molecular Spinorbitals

To perform most of the post-Hartree--Fock calculations, we usually need to transform the integrals to the Molecular Spinorbital basis. We will describe it here and refer to it in the post-Hartree--Fock methods sections. We will also present the post-Hartree--Fock methods using the integrals in the Molecular Spinorbital basis (and its antisymmetrized form in case of the Coulomb integrals), since it is more general.

All the integrals defined in the equations \ref{eq:hamiltonian}, \ref{eq:overlap}, and \ref{eq:coulomb} and even the Fock matrix in the equation \ref{eq:fock} are defined in the basis of atomic orbitals. To transform these integrals to the Molecular Spinorbital basis, we need to use the coefficient matrix $\mathbf{C}$ obtained from the solution of the Roothaan equations \ref{eq:roothaan}. The coefficient matrix $\mathbf{C}$, which is obtained from the Restricted Hartree--Fock calculation, is calculated in the spatial molecular orbital basis. The first step is to expand the coefficient matrix $\mathbf{C}$ to the Molecular Spinorbital basis. This can be done mathematically using the tiling matrix $\mathbf{P}\_{n\times 2n}$, defined as

\begin{equation}
\mathbf{P}=
\begin{pmatrix}
e\_1&e\_1&e\_2&e\_2&\dots&e\_n&e\_n
\end{pmatrix}
,
\end{equation}

where $e\_i$ represents the $i$-th column of the identity matrix $\mathbf{I}\_n$ and the matrices $\mathbf{M}\_{n\times 2n}$ and $\mathbf{N}\_{n\times 2n}$ with elements given by

\begin{equation}
M\_{ij}=1-j\bmod 2,N\_{ij}=j \bmod 2.
\end{equation}

The coefficient matrix $\mathbf{C}$ in the MS basis can be then expressed as

\begin{equation}
\mathbf{C}^{MS}=
\begin{pmatrix}
\mathbf{CP} \\\\\
\mathbf{CP}
\end{pmatrix}
\odot
\begin{pmatrix}
\mathbf{M} \\\\\
\mathbf{N}
\end{pmatrix}
,
\end{equation}

where $\odot$ denotes the Hadamard product. This transformed matrix $\mathbf{C}^{MS}$ is then used to transform the Coulomb integrals $\mathbf{J}$ to the MS basis as

\begin{equation}
J\_{pqrs}^{MS}=C\_{\mu p}^{MS}C\_{\nu q}^{MS}(\mathbf{I}\_{2}\otimes\_K(\mathbf{I}\_{2}\otimes\_K\mathbf{J})^{(4,3,2,1)})\_{\mu\nu\kappa\lambda}C\_{\kappa r}^{MS}C\_{\lambda s}^{MS},
\end{equation}

where the superscript $(4,3,2,1)$ denotes the axes transposition and $\otimes\_K$ is the Kronecker product. This notation accounts for the spin modifications and ensures that the transformations adhere to quantum mechanical principles. We also define the antisymmetrized Coulomb integrals in physicists' notation as

\begin{equation}
\braket{pq||rs}=(J\_{pqrs}^{MS}-J\_{psrq}^{MS})^{(1,3,2,4)}.
\end{equation}

For the transformation of the one-electron integrals such as the core Hamiltonian, the overlap matrix and also the Fock matrix, we use the formula

\begin{equation}
A\_{pq}^{MS}=C\_{\mu p}^{MS}(\mathbf{I}\_{2}\otimes\_K\mathbf{A})\_{\mu\nu}C\_{\nu q}^{MS},
\end{equation}

where $\mathbf{A}$ is an arbitrary one-electron integral.

{:.cite}
> Gill, Peter M. W. 1994. *Molecular Integrals over Gaussian Basis Functions*. Academic Press. <https://doi.org/10.1016/S0065-3276(08)60019-2>.
>

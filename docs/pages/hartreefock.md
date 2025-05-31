---
title: Hartree–Fock Method
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Hartree--Fock Method<!--\label{sec:hartree_fock_method}-->

The Hartree--Fock method is a foundational approach in quantum chemistry, aimed at solving the electronic structure problem in molecules by determining an optimal wavefunction. This method simplifies the complex interactions of electrons through a mean-field approximation, where each electron moves in an average field created by all others. This allows for the use of a single set of orbitals, leading to the construction of the Fock operator and iterative solutions to one-electron equations.

However, Hartree--Fock has notable limitations. Its reliance on a single-determinant wavefunction means it struggles to account for electron correlation. This leads to inaccuracies in energy predictions, particularly for systems with strong electron interactions, such as transition metal complexes or molecules with delocalized electrons.

## Theoretical Background

Our primary objective is to solve the Schrödinger equation in the form

$$
\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi}
\end{equation}
$$

where $$\hat{\mathbf{H}}$$ denotes the molecular Hamiltonian operator, $$\ket{\Psi}$$ represents the molecular wavefunction, and $$E$$ is the total energy of the system. The Hartree--Fock approximates the total wavefunction $$\ket{\Psi}$$ as a single Slater determinant, expressed as

$$
\begin{equation}
\ket{\Psi}=\ket{\chi_1\chi_2\cdots\chi_N}
\end{equation}
$$

where $$\chi_i$$ denotes a spin orbital, and $$N$$ is the total number of electrons. The goal of the Hartree--Fock method is to optimize these orbitals in order to minimize the system's total energy, thereby providing a reliable estimate of the electronic structure.

In the Restricted Hartree--Fock method, we impose a constraint on electron spin, which allows us to work with spatial orbitals instead of spin orbitals. This reformulation expresses the Slater determinant in terms of spatial orbitals as

$$
\begin{equation}
\ket{\Psi}=\ket{\Phi_1\Phi_2\cdots\Phi_{N/2}}
\end{equation}
$$

where $$\Phi_i$$ represents a spatial orbital. Notably, the Restricted Hartree--Fock method requires an even number of electrons to satisfy spin-pairing. In practice, atomic orbitals (whether spin or spatial) are typically expanded in a set of basis functions $$\left\lbrace\phi_i\right\rbrace$$, which are often Gaussian functions, allowing for convenient computation with expansion coefficients. After performing some algebraic manipulations, the Hartree--Fock method leads to the Roothaan equations, which are expressed as

$$
\begin{equation}\label{eq:roothaan}
\mathbf{FC}=\mathbf{SC}\mathbf{\varepsilon}
\end{equation}
$$

where $$\mathbf{F}$$ is the Fock matrix, $$\mathbf{C}$$ is the matrix of orbital coefficients, $$\mathbf{S}$$ is the overlap matrix, and $$\mathbf{\varepsilon}$$ represents the orbital energies. These matrices will be defined in detail later.

## Implementation of the Restricted Hartree--Fock Method

To begin, we define the core Hamiltonian, also known as the one-electron Hamiltonian. This component of the full Hamiltonian excludes electron-electron repulsion and is expressed in index notation as

$$
\begin{equation}\label{eq:hamiltonian}
H_{\mu\nu}^{\mathrm{core}}=T_{\mu\nu}+V_{\mu\nu}
\end{equation}
$$

where $$\mu$$ and $$\nu$$ are indices of the basis functions. Here, $$T_{\mu\nu}$$ represents a kinetic energy matrix element, while $$V_{\mu\nu}$$ denotes a potential energy matrix element. These matrix elements are defined as

$$
\begin{align}
T_{\mu\nu}&=\braket{\phi_{\mu}\vert\hat{T}\vert\phi_{\nu}} \\
V_{\mu\nu}&=\braket{\phi_{\mu}\vert\hat{V}\vert\phi_{\nu}}
\end{align}
$$

Additionally, the overlap integrals, which describe the extent of overlap between basis functions, are given by

$$
\begin{equation}\label{eq:overlap}
S_{\mu\nu}=\braket{\phi_{\mu}\vert\phi_{\nu}}
\end{equation}
$$

and another essential component are the two-electron repulsion integrals, defined as

$$
\begin{equation}\label{eq:coulomb}
J_{\mu\nu\kappa\lambda}=\braket{\phi_{\mu}\phi_{\mu}\vert\hat{J}\vert\phi_{\kappa}\phi_{\lambda}}
\end{equation}
$$

which play crucial roles in the Hartree--Fock calculation. All of these integrals over (Gausssian) basis functions are usually calculated using analytical expressions.<!--\supercite{10.1016/S0065-3276!08!60019-2}--> The solution to the Roothaan equations \eqref{eq:roothaan} requires an iterative procedure, since the Fock matrix defined as

$$
\begin{equation}\label{eq:fock}
F_{\mu\nu}=H_{\mu\nu}^{\mathrm{core}}+D_{\kappa\lambda}\left(J_{\mu\nu\kappa\lambda}-\frac{1}{2}J_{\mu\lambda\kappa\nu}\right)
\end{equation}
$$

depends on the unknown density matrix $$\mathbf{D}$$, defined as

$$
\begin{equation}\label{eq:hf_density}
D_{\mu\nu}=2C_{\mu i}C_{\nu i}
\end{equation}
$$

This iterative procedure is executed through the Self-Consistent Field method. An initial guess is made for the density matrix $$\mathbf{D}$$ (usually zero), the Roothaan equations \eqref{eq:roothaan} are solved and the density matix is updated using the equation \eqref{eq:hf_density}. The total energy of the system is then calculated using the core Hamiltonian and the Fock matrix as

$$
\begin{equation}
E=\frac{1}{2}D_{\mu\nu}(H_{\mu\nu}^{\mathrm{core}}+F_{\mu\nu})+E_{\mathrm{nuc}}
\end{equation}
$$

After convergence of both the density matrix and total energy, the process concludes, yielding the optimized molecular orbitals. The total energy of the system also includes the nuclear repulsion energy, which is given by

$$
\begin{equation}
E_{\mathrm{nuc}}=\sum_{A}\sum_{B<A}\frac{Z_{A}Z_{B}}{R_{AB}}
\end{equation}
$$

where $$Z_A$$ is the nuclear charge of atom $$A$$, and $$R_{AB}$$ is the distance between atoms $$A$$ and $$B$$. The Roothaan equations \eqref{eq:roothaan} are a generalized eigenvalue problem, which can be transformed into a standard eigenvalue problem as explained [here](generalizedeigenvalueproblem.html#generalized-eigenvalue-problem)<!--in Section \ref{sec:generalized_eigenvalue_problem}-->.

### Direct Inversion in the Iterative Subspace

In the Hartree--Fock method, convergence of the density matrix and energy can be significantly accelerated by employing the Direct Inversion in the Iterative Subspace technique. Direct Inversion in the Iterative Subspace achieves this by storing Fock matrices from previous iterations and constructing an optimized linear combination that minimizes the current iteration's error. This approach is especially valuable in Hartree--Fock calculations for large systems, where convergence issues are more common and challenging to resolve.

We start by defining the error vector of $$i$$-th iteration $$\mathbf{e}_i$$ as

$$
\begin{equation}
\mathbf{e}_i=\mathbf{S}_i\mathbf{D}_i\mathbf{F}_i-\mathbf{F}_i\mathbf{D}_i\mathbf{S}_i
\end{equation}
$$

Our goal is to transform the Fock matrix $$\mathbf{F}_i$$ as

$$
\begin{equation}\label{eq:fock_extrapolate}
\mathbf{F}_i=\sum_{j=i-(L+1)}^{i-1}c_j\mathbf{F}_j
\end{equation}
$$

where $$c_j$$ are the coefficients that minimize the error matrix $$\mathbf{e}_i$$ and $$L$$ is the number of Fock matrices we store called the subspace size. To calculate the coefficients $$c_j$$, we solve the set linear equations

$$
\begin{equation}
\begin{bmatrix}
\mathbf{e}_1\cdot\mathbf{e}_1 & \dots & \mathbf{e}_1\cdot\mathbf{e}_{L+1} & 1 \\
\vdots & \ddots & \vdots & \vdots \\
\mathbf{e}_{L+1}\cdot\mathbf{e}_1 & \dots & \mathbf{e}_{L+1}\cdot\mathbf{e}_{L+1} & 1 \\
1 & \dots & 1 & 0 \\
\end{bmatrix}
\begin{bmatrix}
c_1 \\
\vdots \\
c_{L+1} \\
\lambda \\
\end{bmatrix}
=
\begin{bmatrix}
0 \\
\vdots \\
0 \\
1 \\
\end{bmatrix}
\end{equation}
$$

After solving the linear equations, we use the coefficients $$c_j$$ to construct the new Fock matrix $$\mathbf{F}_i$$ according to the equation \eqref{eq:fock_extrapolate} and proceed as usual with the Hartree--Fock calculation.

### Gradient of the Restricted Hartree--Fock Method

If we perform the calculation as described above and get the density matrix $$\mathbf{D}$$ we can evaluate the nuclear energy gradient as<!--\supercite{10.1002/9780470749593.hrs006}-->

$$
\begin{equation}\label{eq:hf_gradient}
\frac{\partial E}{\partial X_{A,i}}=\sum_{\mu\nu\in\left\lbrace\phi_A\right\rbrace}D_{\mu\nu}\frac{\partial H_{\mu\nu}^{\mathrm{core}}}{\partial X_{A,i}}+2\sum_{\mu\nu\in\left\lbrace\phi_A\right\rbrace}D_{\mu\nu}D_{\kappa\lambda}\frac{\partial\left(J_{\mu\nu\kappa\lambda}-\frac{1}{2}J_{\mu\lambda\kappa\nu}\right)}{\partial X_{A,i}}-2W_{\mu\nu}\frac{\partial S_{\mu\nu}}{\partial X_{A,i}}
\end{equation}
$$

where $$A$$ is the index of an atom, $$i$$ is the index of the coordinate, $$\left\lbrace\phi_A\right\rbrace$$ is the set of all basis functions located at atom $$A$$ and $$\mathbf{W}$$ is energy weighed density matrix defined as

$$
\begin{equation}
W_{\mu\nu}=2C_{\mu i}C_{\nu i}\varepsilon_i
\end{equation}
$$

Keep in mind that the indices $$\kappa$$ and $$\lambda$$ in the gradient equation \eqref{eq:hf_gradient} are summed over all basis functions.

## Integral Transforms to the Basis of Molecular Spinorbitals<!--\label{sec:integral_transform}-->

To carry out most post-Hartree--Fock calculations, it is essential to transform the integrals into the Molecular Spinorbital basis. We will outline this transformation process here and refer to it in subsequent sections on post-Hartree--Fock methods. The post-Hartree--Fock methods in this document will be presented using the integrals in the Molecular Spinorbital basis (and in its antisymmetrized form for the two-electron integrals) as this approach is more general.

All the integrals defined in the equations \eqref{eq:hamiltonian}, \eqref{eq:overlap}, and \eqref{eq:coulomb}, as well as Fock matrix in the equation \eqref{eq:fock} are defined in the basis of atomic orbitals. To transform these integrals to the Molecular Spinorbital basis, we utilize the coefficient matrix $$\mathbf{C}$$ obtained from the solution of the Roothaan equations \eqref{eq:roothaan}. This coefficient matrix $$\mathbf{C}$$ is initially calculated in the spatial molecular orbital basis (in the Restricted Hartree--Fock calculation).

The first step involves expanding the coefficient matrix $$\mathbf{C}$$ to the Molecular Spinorbital basis. This transformation can be mathematically expressed using the tiling matrix $$\mathbf{P}_{n\times 2n}$$, defined as

$$
\begin{equation}
\mathbf{P}=
\begin{pmatrix}
e_1&e_1&e_2&e_2&\dots&e_n&e_n
\end{pmatrix}
\end{equation}
$$

where $$e_i$$ represents the $$i$$-th column of the identity matrix $$\mathbf{I}_n$$. Additionally, we define the matrices $$\mathbf{M}_{n\times 2n}$$ and $$\mathbf{N}_{n\times 2n}$ with elements given by

$$
\begin{equation}
M_{ij}=1-j\bmod 2,N_{ij}=j \bmod 2
\end{equation}
$$

The coefficient matrix $$\mathbf{C}$$ in the Molecular Spinorbital basis can be then expressed as

$$
\begin{equation}
\mathbf{C}^{\mathrm{MS}}=
\begin{pmatrix}
\mathbf{CP} \\
\mathbf{CP}
\end{pmatrix}
\odot
\begin{pmatrix}
\mathbf{M} \\
\mathbf{N}
\end{pmatrix}
\end{equation}
$$

where $$\odot$$ denotes the Hadamard product. This transformed matrix $$\mathbf{C}^{\mathrm{MS}}$$ is subsequently used to transform the two-electron integrals $$\mathbf{J}$$ to the Molecular Spinorbital basis as

$$
\begin{equation}
J_{pqrs}^{\mathrm{MS}}=C_{\mu p}^{\mathrm{MS}}C_{\nu q}^{\mathrm{MS}}(\mathbf{I}_{2}\otimes_K(\mathbf{I}_{2}\otimes_K\mathbf{J})^{(4,3,2,1)})_{\mu\nu\kappa\lambda}C_{\kappa r}^{\mathrm{MS}}C_{\lambda s}^{\mathrm{MS}}
\end{equation}
$$

where the superscript $$(4,3,2,1)$$ denotes the axes transposition and $$\otimes_K$$ is the Kronecker product. This notation accommodates the spin modifications and ensures adherence to quantum mechanical principles. We also define the antisymmetrized two-electron integrals in physicists' notation as

$$
\begin{equation}
\braket{pq\vert\vert rs}=(J_{pqrs}^{\mathrm{MS}}-J_{psrq}^{\mathrm{MS}})^{(1,3,2,4)}
\end{equation}
$$

For the transformation of the one-electron integrals such as the core Hamiltonian, the overlap matrix and also the Fock matrix, we use the formula

$$
\begin{equation}
A_{pq}^{\mathrm{MS}}=C_{\mu p}^{\mathrm{MS}}(\mathbf{I}_{2}\otimes_K\mathbf{A})_{\mu\nu}C_{\nu q}^{\mathrm{MS}}
\end{equation}
$$

where $$\mathbf{A}$$ is an arbitrary matrix of one-electron integrals. Since many post-Hartree--Fock methods rely on differences of orbital energies in the denominator, we define the tensors

$$
\begin{align}
\varepsilon^{a}_{i}&=\varepsilon_i-\varepsilon_a \\
\varepsilon^{ab}_{ij}&=\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b \\
\varepsilon^{abc}_{ijk}&=\varepsilon_i+\varepsilon_j+\varepsilon_k-\varepsilon_a-\varepsilon_b-\varepsilon_c
\end{align}
$$

for convenience. These tensors enhance code readability and efficiency, making it easier to understand and work with the underlying mathematical framework. Here and also throughout the rest of the document, the indices $$i$$, $$j$$ and $$k$$ run over occupied orbitals, whereas the indices $$a$$, $$b$$ and $$c$$ run over virtual orbitals.

{:.cite}
> Gill, Peter M. W. 1994. *Molecular Integrals over Gaussian Basis Functions*. Academic Press. <https://doi.org/10.1016/S0065-3276(08)60019-2>.
>

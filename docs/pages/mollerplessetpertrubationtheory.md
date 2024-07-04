---
title: Møller–Plesset Pertrubation Theory
parent: Electronic Structure Methods
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Møller–Plesset Perturbation Theory (MPPT)

Møller-Plesset perturbation theory is a quantum mechanical method used to improve the accuracy of electronic structure calculations within the framework of Hartree-Fock theory. It involves treating electron-electron correlation effects as a perturbation to the reference Hartree-Fock wave function. The method is named after its developers, physicists C. Møller and M. S. Plesset. By systematically including higher-order corrections, Møller-Plesset perturbation theory provides more accurate predictions of molecular properties compared to the initial Hartree-Fock approximation.

## Møller–Plesset Perturbation Theory of n-th Order

Møller–Plesset perturbation theory of n-th order (MPN) specifically considers electron correlation effects by introducing perturbative corrections up to the n-th order.

### The Transform of Integrals to Molecular Spinorbital Basis

To begin, we need to convert the coefficient matrix $\mathbf{C}$ into the basis of molecular spinorbitals (MS), as this matrix is essential for later transformations of the Coulomb integrals $\mathbf{J}$. Let the number of spatial basis functions be denoted as $n$. The coefficient matrix in the MS basis will then be of dimension $2n \times 2n$, to incorporate the electron spin. In this transformation, each column (orbital) is duplicated and the result is vertically concatenated with itself. The corresponding segments are then zeroed out to respect the Pauli exclusion principle. This is done mathematically using the tiling matrix $\mathbf{P}_{n \times 2n}$, defined as

\begin{equation}
\mathbf{P}=\begin{pmatrix}e_1&e_1&e_2&e_2&\dots&e_n&e_n\end{pmatrix}
\end{equation}

where $e_i$ represents the $i$-th column of the identity matrix $\mathbf{I}\_n$. To facilitate the correct allocation of spin states, we define the matrices $\mathbf{M}\_{n \times 2n}$ and $\mathbf{N}_{n \times 2n}$ with elements given by

\begin{equation}
M_{ij}=1-j\bmod 2,N_{ij}=j \bmod 2
\end{equation}

The coefficient matrix $\mathbf{C}$ in the MS basis is then expressed as

\begin{equation}
\mathbf{C}^{MS}=\begin{pmatrix}\mathbf{CP}\\\ \mathbf{CP}\end{pmatrix}\odot\begin{pmatrix}\mathbf{M}\\\ \mathbf{N}\end{pmatrix}
\end{equation}

where $\odot$ denotes the Hadamard product. This transformed matrix $\mathbf{C}^{MS}$ is then used to transform the Coulomb integrals to the MS basis. The Coulomb integrals in the MS basis are given by

\begin{equation}
J_{pqrs}^{MS}=C_{\mu p}^{MS}C_{\nu q}^{MS}(\mathbf{I}\_{2}\otimes_K(\mathbf{I}\_{2}\otimes_K\mathbf{J})^{(4,3,2,1)})\_{\mu\nu\kappa\lambda}C_{\kappa r}^{MS}C_{\lambda s}^{MS}
\end{equation}

where the superscript $(4,3,2,1)$ denotes the axes transposition. This notation accounts for the spin modifications and ensures that the transformations adhere to quantum mechanical principles. Let's also define the antisymmetrized Coulomb integrals in the MS basis as

\begin{equation}
[pq||rs]=J_{pqrs}^{MS}-J_{psrq}^{MS}
\end{equation}

These integrals are calculated just to simplify the notation in the following formulas.

### The Double Excitation Amplitudes

For the definition of general n-th order formulas we well need an expression for the double excitation amplitudes, which are defined as

\begin{equation}
t_{ij}^{ab}=\frac{[ia||jb]}{\varepsilon_a+\varepsilon_b-\varepsilon_i-\varepsilon_j}
\end{equation}

where $a$ and $b$ are unoccupied (virtual) spinorbitals and $i$ and $j$ are occupied spinorbitals. The vector $\varepsilon$ encompasses energies of spinorbitals, obtained through the solution of the Roothaan equations.

### Correlation Energy Calculation

Having computed the Coulomb integral in the molecular spinorbital (MS) basis and the double excitation amplitudes, we can now calculate the correlation energy for the Møller–Plesset perturbation theory.

The 2nd-order correlation energy:

\begin{equation}
E_{corr}^{MP2}=-\frac{1}{4}\sum_{ijab}t_{ij}^{ab}[ia||jb]
\end{equation}

The 3rd-order correlation energy:

\begin{equation}
E_{corr}^{MP3}=\frac{1}{8}\sum_{ijab}t_{ij}^{ab}[ac||bd]t_{ij}^{cd}+\frac{1}{8}\sum_{ijab}t_{ij}^{ab}[ki||lj]t_{kl}^{ab}-\sum_{ijab}t_{ij}^{ab}[ki||bc]t_{kj}^{ac}
\end{equation}

Remember that the sums over $i$, $j$, $k$, and $l$ are over occupied spinorbitals, while the sums over $a$, $b$, $c$, and $d$ are over virtual spinorbitals.

{:.note}
> Keep in mind that for the MP2 method we could integrate out the spin and obtain the correlation energy as
> \begin{equation}
> E_{corr}^{MP2}=\sum_{iajb}\frac{J_{iajb}^{MO}(2J_{iajb}^{MO}-J_{ibja}^{MO})}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
> \end{equation}
> where $J_{pqrs}^{MO}$ are the Coulomb integrals in the molecular orbital (MO) basis calculated as
> \begin{equation}
> J_{pqrs}^{MO}=C_{\mu p}C_{\nu q}J_{\mu\nu\kappa\lambda}C_{\kappa r}C_{\lambda s}
> \end{equation}
> and $\varepsilon_i$ are the spatial orbital energies obtained from the solution of the Roothaan equations.

---
title: Møller–Plesset Pertrubation Theory
parent: Electronic Structure Methods
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Møller–Plesset Perturbation Theory

Møller-Plesset perturbation theory is a quantum mechanical method used to improve the accuracy of electronic structure calculations within the framework of Hartree-Fock theory. It involves treating electron-electron correlation effects as a perturbation to the reference Hartree-Fock wave function. The method is named after its developers, physicists C. Møller and M. S. Plesset. By systematically including higher-order corrections, Møller-Plesset perturbation theory provides more accurate predictions of molecular properties compared to the initial Hartree-Fock approximation.

## Restricted Møller–Plesset Perturbation Theory of 2nd Order

Restricted Møller–Plesset perturbation theory of 2nd order (RMP2) specifically considers electron correlation effects by introducing perturbative corrections up to the second order. This method is applicable to closed-shell systems, where the spatial orbitals are occupied by a fixed number of electrons, and it provides a computationally efficient way to account for electron correlation and improve the accuracy of predicted molecular properties.

### Coulomb Integral Transform

Prior to computing the correlation energy, it is necessary to convert the Coulomb integral $\mathbf{J}$ into the basis of molecular orbitals (MO). This transformation is straightforward and is executed using the formula

\begin{equation}
J_{pqrs}^{MO}=C_{\mu p}C_{\nu q}J_{\mu\nu\kappa\lambda}C_{\kappa r}C_{\lambda s}
\end{equation}

where $\mathbf{C}$ represents the matrix of coefficients obtained from the Hartree–Fock method.

### Correlation Energy Calculation

Having computed the Coulomb integral in the molecular orbital (MO) basis, the correlation energy of a closed-shell system can be determined by assessing the sum

\begin{equation}
E_{corr}^{MP2}=\sum_{iajb}\frac{J_{iajb}^{MO}(2J_{iajb}^{MO}-J_{ibja}^{MO})}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
\end{equation}

where indices $i$ and $j$ refer to occupied spatial orbitals, and $a$ and $b$ refer to unoccupied spatial orbitals. The vector $\varepsilon$ encompasses orbital energies, obtained through the solution of the Roothaan equations.

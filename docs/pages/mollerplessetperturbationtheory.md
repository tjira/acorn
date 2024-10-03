---
title: Møller–Plesset Perturbation Theory
parent: Electronic Structure Methods
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Møller--Plesset Perturbation Theory

Møller--Plesset Perturbation Theory is a quantum mechanical method used to improve the accuracy of electronic structure calculations within the framework of Hartree--Fock theory. It involves treating electron-electron correlation effects as a perturbation to the reference Hartree--Fock wave function. The method is named after its developers, physicists C. Møller and M. S. Plesset. By systematically including higher-order corrections, Møller--Plesset Perturbation Theory provides more accurate predictions of molecular properties compared to the initial Hartree--Fock approximation.

## Theory of the Pertrubative Approach

As for the Hartree--Fock method, we start with the Schrödinger equation in the form

\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi},
\end{equation}

where $\hat{\mathbf{H}}$ is the molecular Hamiltonian operator, $\ket{\Psi}$ is the molecular wave function, and $E$ is the total energy of the system. In the Møller--Plesset perturbation theory we write the Hamiltonian operator as

\begin{equation}
\hat{\mathbf{H}}=\hat{\mathbf{H}}^{(0)}+\lambda\hat{\mathbf{H}}^{'},
\end{equation}

where $\hat{\mathbf{H}}^{(0)}$ is the Hamiltonian used in the Hartree--Fock method (the electrons are moving in the mean field), $\lambda$ is a parameter between 0 and 1, and $\hat{\mathbf{H}}^{'}$ is the perturbation operator representing the missing electron-electron interactions. We can then expand the wavefunction $\ket{\Psi}$ and total energy $E$ in a power series of $\lambda$ as

\begin{align}
\ket{\Psi}&=\ket{\Psi^{(0)}}+\lambda\ket{\Psi^{(1)}}+\lambda^2\ket{\Psi^{(2)}}+\dots \\\\\
E&=E^{(0)}+\lambda E^{(1)}+\lambda^2 E^{(2)}+\dots
\end{align}

and ask, how how does the total energy change with the included terms. After some algebra, we can show that the first-order correction to the total energy is zero, the second-order correction is given by

\begin{equation}
E_{corr}^{MP2}=\sum\_{s>0}\frac{H\_{0s}^{'}H\_{s0}^{'}}{E\_0-E\_s},
\end{equation}

where $s$ runs over all doubly excited determinants, $H\_{0s}^{'}$ is the matrix element of the perturbation operator between the HF determinant and the doubly excited determinant, and $E\_0$ and $E\_s$ are the energies of the reference and doubly excited determinants, respectively.<!--\cite{10.1002/wcms.58}--> We could express all higher-order corrections in a similar way, using only the matrix elements of the perturbation operator and the energies of the determinants. For practical calculations, we apply Slater-Condon rules to evaluate the matrix elements and use the orbital energies obtained from the Hartree--Fock calculation. The expressions for calculation are summarised below.

## Implementation of 2nd and 3rd Order Corrections

Having the antisymmetrized Coulomb integrals in the MS basis and physicists' notation defined [here](hartreefockmethod.html#the-integral-transforms), we can now proceed with the calculation of the correlation energy. We wil use the convention, that the indices $i$, $j$, $k$, and $l$ run over occupied spinorbitals, while the indices $a$, $b$, $c$, and $d$ run over virtual spinorbitals. The 2nd-order and 3rd-order correlation energies are given by the following expressions.

The 2nd-order correlation energy:<!--\cite{10.1021/acs.jpclett.9b01641}-->

\begin{equation}
E_{corr}^{MP2}=\frac{1}{4}\sum_{ijab}\frac{\braket{ab||ij}\braket{ij||ab}}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
\end{equation}

The 3rd-order correlation energy:<!--\cite{10.1021/acs.jpclett.9b01641}-->

\begin{align}
E_{corr}^{MP3}=&\frac{1}{8}\sum_{ijab}\frac{\braket{ab||ij}\braket{cd||ab}\braket{ij||cd}}{(\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b)(\varepsilon_i+\varepsilon_j-\varepsilon_c-\varepsilon_d)}+\nonumber \\\\\
&+\frac{1}{8}\sum_{ijab}\frac{\braket{ab||ij}\braket{ij||kl}\braket{kl||ab}}{(\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b)(\varepsilon_k+\varepsilon_l-\varepsilon_a-\varepsilon_b)}+\nonumber \\\\\
&+\sum_{ijab}\frac{\braket{ab||ij}\braket{cj||kb}\braket{ik||ac}}{(\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b)(\varepsilon_i+\varepsilon_k-\varepsilon_a-\varepsilon_c)}
\end{align}

To calculate the 4th order correction, we would need to write 39 terms, which is not practical. Higher-order corrections are usually not programmed this way, instead, the diagrammatic approach is used.

{:.note}
> Keep in mind that for the MP2 method we could integrate out the spin and obtain the correlation energy as
> \begin{equation}
> E_{corr}^{MP2}=\frac{1}{4}\sum_{iajb}\frac{J_{iajb}^{MO}(2J_{iajb}^{MO}-J_{ibja}^{MO})}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
> \end{equation}
> where $J_{pqrs}^{MO}$ are the Coulomb integrals in the molecular orbital (MO) basis calculated as
> \begin{equation}
> J_{pqrs}^{MO}=C_{\mu p}C_{\nu q}J_{\mu\nu\kappa\lambda}C_{\kappa r}C_{\lambda s}
> \end{equation}
> and $\varepsilon_i$ are the spatial orbital energies obtained from the solution of the Roothaan equations.

{:.cite}
> Cremer, Dieter. 2011. “Møller–Plesset Perturbation Theory: From Small Molecule Methods to Methods for Thousands of Atoms.” *WIREs Computational Molecular Science* 1: 509–30. <https://doi.org/10.1002/wcms.58>.
>
> Bertels, Luke W., Joonho Lee, and Martin Head-Gordon. 2019. “Third-Order Møller–Plesset Perturbation Theory Made Useful? Choice of Orbitals and Scaling Greatly Improves Accuracy for Thermochemistry, Kinetics, and Intermolecular Interactions.” *The Journal of Physical Chemistry Letters* 10: 4170–76. <https://doi.org/10.1021/acs.jpclett.9b01641>.
>

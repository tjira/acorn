---
title: Møller–Plesset Perturbation Theory
parent: Electronic Structure Methods
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Møller--Plesset Perturbation Theory<!--\label{sec:moller_plesset_perturbation_theory}-->

Møller--Plesset Perturbation Theory is a quantum mechanical method used to improve the accuracy of electronic structure calculations within the framework of Hartree--Fock theory. It involves treating electron-electron correlation effects as a perturbation to the reference Hartree--Fock wave function. The method is named after its developers, physicists C. Møller and M. S. Plesset. By systematically including higher-order corrections, Møller--Plesset Perturbation Theory provides more accurate predictions of molecular properties compared to the initial Hartree--Fock approximation.

## Theory of the Perturbative Approach

As for the Hartree--Fock method, we start with the Schrödinger equation in the form

\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi}
\end{equation}

where $\hat{\mathbf{H}}$ is the molecular Hamiltonian operator, $\ket{\Psi}$ is the molecular wave function, and $E$ is the total energy of the system. In the Møller--Plesset perturbation theory we write the Hamiltonian operator as

\begin{equation}
\hat{\mathbf{H}}=\hat{\mathbf{H}}^{(0)}+\lambda\hat{\mathbf{H}}^{'}
\end{equation}

where $\hat{\mathbf{H}}^{(0)}$ is the Hamiltonian used in the Hartree--Fock method (representing electrons moving in the mean field), $\lambda$ is a parameter between 0 and 1, and $\hat{\mathbf{H}}^{'}$ is the perturbation operator representing the missing electron-electron interactions not included in the Hartree--Fock approximation. We then expand the wavefunction $\ket{\Psi}$ and total energy $E$ as a power series in $\lambda$ as

\begin{align}
\ket{\Psi}&=\ket{\Psi^{(0)}}+\lambda\ket{\Psi^{(1)}}+\lambda^2\ket{\Psi^{(2)}}+\dots \\\\\
E&=E^{(0)}+\lambda E^{(1)}+\lambda^2 E^{(2)}+\dots
\end{align}

and ask, how how does the total energy change with the included terms. After some algebra, we can show that the first order correction to the total energy is zero, the second order correction is given by

\begin{equation}
E_{\mathrm{corr}}^{\mathrm{MP2}}=\sum\_{s>0}\frac{H\_{0s}^{'}H\_{s0}^{'}}{E\_0-E\_s}
\end{equation}

where $s$ runs over all doubly excited determinants, $H\_{0s}^{'}$ is the matrix element of the perturbation operator between the Hartree--Fock determinant and the doubly excited determinant, and $E\_0$ and $E\_s$ are the energies of the reference and doubly excited determinants, respectively.<!--\supercite{10.1002/wcms.58,1014569052}--> We could express all higher-order corrections in a similar way, using only the matrix elements of the perturbation operator and the energies of the determinants. For practical calculations, we apply Slater-Condon rules to evaluate the matrix elements and use the orbital energies obtained from the Hartree--Fock calculation. The expressions for calculation are summarised below.

## Implementation of 2nd and 3rd Order Corrections

Having the antisymmetrized two-electron integrals in the Molecular Spinorbital basis and physicists' notation defined [here](hartreefock.html#integral-transforms-to-the-basis-of-molecular-spinorbitals)<!--in Section \ref{sec:integral_transform}-->, we can now proceed with the calculation of the correlation energy. The 2nd order correlation energy can be expressed as

\begin{equation}
E_{\mathrm{corr}}^{\mathrm{MP2}}=\frac{1}{4}\sum\_{ijab}\frac{\braket{ab||ij}\braket{ij||ab}}{\varepsilon\_{ij}^{ab}}
\end{equation}

and the 3rd order correlation energy as

\begin{align}
E_{\mathrm{corr}}^{\mathrm{MP3}}=&\frac{1}{8}\sum\_{ijab}\frac{\braket{ab||ij}\braket{cd||ab}\braket{ij||cd}}{\varepsilon\_{ij}^{ab}\varepsilon\_{ij}^{cd}}+\nonumber \\\\\
&+\frac{1}{8}\sum_{ijab}\frac{\braket{ab||ij}\braket{ij||kl}\braket{kl||ab}}{\varepsilon\_{ij}^{ab}\varepsilon\_{kl}^{ab}}+\nonumber \\\\\
&+\sum_{ijab}\frac{\braket{ab||ij}\braket{cj||kb}\braket{ik||ac}}{\varepsilon\_{ij}^{ab}\varepsilon\_{ik}^{ac}}
\end{align}

To calculate the 4th order correction, we would need to write 39 terms, which is not practical. Higher-order corrections are usually not programmed this way, instead, the diagrammatic approach is used.<!--\supercite{1014569052,10.1016/0010-4655!73!90016-7,10.1016/0010-4655!73!90017-9}-->

{:.cite}
> Cremer, Dieter. 2011. “Møller–Plesset Perturbation Theory: From Small Molecule Methods to Methods for Thousands of Atoms.” *WIREs Computational Molecular Science* 1: 509–30. <https://doi.org/10.1002/wcms.58>.
>
> Paldus, Josef, and Heymans H. C. Wong. 1973. “Computer Generation of Feynman Diagrams for Perturbation Theory i. General Algorithm.” *Computer Physics Communications* 6: 1–7. <https://doi.org/10.1016/0010-4655(73)90016-7>.
>
> Wong, Heymans H. C., and Josef Paldus. 1973. “Computer Generation of Feynman Diagrams for Perturbation Theory II. Program Description.” *Computer Physics Communications* 6: 9–16. <https://doi.org/10.1016/0010-4655(73)90017-9>.
>
> Szabo, Attila, and Neil S. Ostlund. 1996. *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*. Courier Corporation. <https://www.worldcat.org/title/1014569052>.
>

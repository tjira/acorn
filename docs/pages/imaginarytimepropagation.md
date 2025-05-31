---
title: Imaginary Time Propagation
parent: Time Evolution in Quantum Mechanics
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Imaginary Time Propagation

Imaginary time propagation is a numerical technique used to find the ground state of a quantum system by evolving it in imaginary time instead of real time. By replacing time with an imaginary variable, higher energy states decay exponentially, leaving only the lowest-energy state. This method is widely used in quantum chemistry and condensed matter physics to determine electronic and nuclear ground states. However, like real time propagation, its applicability is limited to small systems due to the high computational cost of solving the underlying equations.

## Theoretical Background

Similar to real time propagation [here](rtp.html), we start with the Time-Dependent Schrödinger Equation, but this time we perform a suspicious substitution $\tau\rightarrow\mathrm{i}t$ to obtain the equation

\begin{equation}\label{eq:tdse_it}
\hbar\frac{\mathrm{d}}{\mathrm{d}\tau}\ket{\Psi\left(\tau\right)}=-\hat{\mathbf{H}}\ket{\Psi\left(\tau\right)},
\end{equation}

which has a formal solution given by

\begin{equation}\label{eq:wf_it}
\ket{\Psi\left(\tau\right)}=\mathrm{e}^{-\frac{\hat{\mathbf{H}}\tau}{\hbar}}\ket{\Psi\left(0\right)}.
\end{equation}

We now expand the wavefunction in the basis of the eigenstates of the Hamiltonian $\hat{\mathbf{H}}$ as

\begin{equation}\label{eq:wf_expansion_eigenstates}
\ket{\Psi\left(\tau\right)}=\sum_{n}c_{n}\left(\tau\right)\ket{\phi_{n}},
\end{equation}

where $\ket{\phi_{n}}$ are the eigenstates of the Hamiltonian $\hat{\mathbf{H}}$ and $c_{n}\left(\tau\right)$ are the expansion coefficients. Substituting the wavefunction expression \eqref{eq:wf_expansion_eigenstates} into the solution of Time-Dependent Schrödinger Equation in imaginary time \eqref{eq:wf_it} we obtain

\begin{equation}\label{eq:tdse_it_expansion}
\ket{\Psi\left(\tau\right)}=\sum_{n}c_{n}\left(\tau\right)\mathrm{e}^{-\frac{E_{n}\tau}{\hbar}}\ket{\phi_{n}},
\end{equation}

where $E_{n}$ are the eigenvalues of the Hamiltonian $\hat{\mathbf{H}}$. We now see, that the bigger the value of $E_n$, the faster the corresponding term in the sum decays. This means that, in the limit $\tau\rightarrow\infty$, only the term corresponding to the smallest eigenvalue $E_{0}$ will survive, and the wavefunction will converge to the ground state of the system. Note that we need to normalize the wavefunction at each step of the imaginary time propagation to ensure that the wavefunction remains normalized, otherwise the wavefunction will decay to zero.

---
title: Real Time Propagation
parent: Quantum Dynamics
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Real Time Propagation

Real time propagation is a computational approach used to study the time evolution of quantum systems according to the Time Dependent Schrödinger Equation. It is essential for modeling dynamic processes such as electronic excitations, molecular vibrations, and non equilibrium phenomena in quantum chemistry. By solving the wavefunction’s evolution step by step in real time, this method captures transient states and ultrafast reactions. However, due to the exponential growth of computational complexity, real time propagation is feasible only for small systems or simplified models.

## Theoretical Background

The time evolution in quantum mechanics is governed by the Time Dependent Schrödinger Equation in the form

\begin{equation}
i\hbar\frac{\mathrm{d}}{\mathrm{d} t}\ket{\Psi\left(t\right)}=\hat{\mathbf{H}}\ket{\Psi\left(t\right)},
\end{equation}

where $i$ is the imaginary unit, $\hbar$ is the reduced Planck constant, $\ket{\Psi\left(t\right)}$ is the time dependent wavefunction, $\hat{\mathbf{H}}$ is the Hamiltonian and  $x$ represents the position of a quantum particle, such as an electron or proton. The solution of the Time Dependent Schrödinger Equation can be in general written as

\begin{equation}
\ket{\Psi\left(t\right)}=\hat{\mathbf{U}}\left(t\right)\ket{\Psi\left(0\right)},
\end{equation}

where the operator $\hat{\mathbf{U}}$ is called a time evolution operator. Application of the propagator to a wavefunction is equivalent to evolving the wave function from time 0 to time $t$. The knowledge of the evolution operator is, therefore, critical for a description of dynamics in quantum systems.

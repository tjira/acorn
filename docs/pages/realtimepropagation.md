---
title: Real Time Propagation
parent: Quantum Dynamics
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Real Time Propagation<!--\label{sec:real_time_propagation}-->

Real time propagation is a computational approach used to study the time evolution of quantum systems according to the Time-Dependent Schrödinger Equation. It is essential for modeling dynamic processes such as electronic excitations, molecular vibrations, and non equilibrium phenomena in quantum chemistry. By solving the wavefunction’s evolution step by step in real time, this method captures transient states and ultrafast reactions. However, due to the exponential growth of computational complexity, real time propagation is feasible only for small systems or simplified models.

## Theoretical Background

The time evolution in quantum mechanics is governed by the Time-Dependent Schrödinger Equation in the form

\begin{equation}\label{eq:tdse}
\mathrm{i}\hbar\frac{\mathrm{d}}{\mathrm{d}t}\ket{\Psi\left(t\right)}=\hat{\mathbf{H}}\ket{\Psi\left(t\right)},
\end{equation}

where $i$ is the imaginary unit, $\hbar$ is the reduced Planck constant, $\ket{\Psi\left(t\right)}$ is the time dependent wavefunction, $\hat{\mathbf{H}}$ is the Hamiltonian operator, which encodes the total energy of the system. In the position representation, one could write $\Psi\left(x,t\right)$, where $x$ is the position coordinate for a quantum particle such as an electron or proton.

A general solution of the Time-Dependent Schrödinger Equation can be written as

\begin{equation}\label{eq:tdse_propagator}
\ket{\Psi\left(t\right)}=\hat{\mathbf{U}}\left(t\right)\ket{\Psi\left(0\right)},
\end{equation}

where the operator $\hat{\mathbf{U}}$ is called a time evolution operator (propagator). The action if $\hat{\mathbf{U}}$ on the initial wavefunction $\ket{\Psi\left(0\right)}$ propagates the state from time $0$ to time $t$, carrying all the information about the time evolution of the system. The propagator is straightforward to find if the Hamiltonian $\hat{\mathbf{H}}$ does not depend on time. In this case, the formal solution to the Time-Dependent Schrödinger Equation can be written as

\begin{equation}\label{eq:tdse_explicit_propagator}
\ket{\Psi\left(t\right)}=\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}t}\ket{\Psi\left(0\right)}.
\end{equation}

Comparing equations \eqref{eq:tdse_propagator} and \eqref{eq:tdse_explicit_propagator} reveals that the time evolution operator for a time independent Hamiltonian can be expressed as

\begin{equation}\label{eq:propagator}
\hat{\mathbf{U}}(t)=\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}t}.
\end{equation}

This exponential operator contains the complete description of how the system evolves with time. In practical applications, however, evaluating $\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}t}$ analytically is rarely straightforward for nontrivial Hamiltonians, which often necessitates numerical methods for an accurate treatment of quantum dynamics.

### Space Discretization and the Diabatic Basis

To represent the wavefunction numerically, we sample it on a spatial grid. In practical terms, this means introducing a dependence on the coordinate $x$, so that $\ket{\Psi\left(t\right)}$ becomes $\ket{\Psi\left(x,t\right)}$. More importantly, we expand the wavefunction $\ket{\Psi\left(x,t\right)}$ in a set of basis functions $\lbrace\Phi_i\rbrace_{i=1}^N$, where $N$ is the number of states included in the dynamics. In the following derivations, we will consider $N=2$, but the extension to multiple states (or even one state) is straightforward. The basis set will be orthonormal and satisfy the condition

\begin{equation}\label{eq:diabatic_condition}
\bra{\Phi_i}\frac{\partial}{\partial x}\ket{\Phi_j}=0,\,i\neq j.
\end{equation}

This condition essentially means that the states $\ket{\Phi_i}$ are either constant in space, or purely imaginary and a basis with this property is called diabatic diabatic basis. Note that the specific functional form of these basis states is not crucial for the dynamics, since they are implicitly taken to be the canonical basis in $\mathbb{R}^N$ when defining the potential energy surfaces. For a two-state system, our wavefunction in the diabatic basis can be written as

\begin{equation}\label{eq:diabatic_wavefunction}
\ket{\Psi\left(x,t\right)}=c_1\left(x,t\right)\ket{\Phi_1}+c_2\left(x,t\right)\ket{\Phi_2}=
\begin{pmatrix}
c_1\left(x,t\right) \\\\\
c_2\left(x,t\right)
\end{pmatrix},
\end{equation}

where $c_i\left(x,t\right)$ are the expansion coefficients of the wavefunction in the diabatic basis. The Hamiltonian operator $\hat{\mathbf{H}}$ in the diabatic basis takes the form

\begin{equation}\label{eq:diabatic_hamiltonian}
\hat{\mathbf{H}}=\hat{\mathbf{T}}+\hat{\mathbf{V}}=-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}\mathbf{I}+
\begin{pmatrix}
V_{11}(x) & V_{12}(x) \\\\\
V_{21}(x) & V_{22}(x)
\end{pmatrix},
\end{equation}

where $\mathbf{I}$ is the identity matrix, $\hat{\mathbf{T}}$ is the kinetic energy operator, and $\mathbf{V}$ is the potential energy matrix. The diagonal elements $V_{ii}(x)$ are the potential energy surfaces for the diabatic state $i$, while the off-diagonal elements are coupling between the two states.

### Wavefunction Propagation

We aim to reduce the action of the propagator \eqref{eq:propagator} on the diabatic wavefunction \eqref{eq:diabatic_wavefunction} to a simple matrix multiplication. The main difficulty is the differential form of the kinetic operator in the equation \ref{eq:diabatic_hamiltonian}. To handle this, we use a Fourier-based method for applying the kinetic operator $\hat{\mathbf{T}}$ on a wavefunction $\ket{\Psi\left(x,t\right)}$. Taking the Fourier Transform $\hat{\mathbf{T}}\ket{\Psi\left(x,t\right)}$ yields

\begin{align}\label{eq:kinetic_fourier_method}
\mathcal{F}\left\lbrace\hat{\mathbf{T}}\ket{\Psi\left(x,t\right)}\right\rbrace&=\mathcal{F}\left\lbrace-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}\ket{\Psi\left(x,t\right)}\right\rbrace=\frac{\hbar^2 k^2}{2m}\mathcal{F}\left\lbrace\ket{\Psi\left(x,t\right)}\right\rbrace \\\\\
\hat{\mathbf{T}}\ket{\Psi\left(x,t\right)}&=\frac{\hbar^2}{2m}\mathcal{F}^{-1}\left\lbrace k^2\mathcal{F}\left\lbrace\ket{\Psi\left(x,t\right)}\right\rbrace\right\rbrace,
\end{align}

where $k$ is the coordinate in the Fourier space. In other words, to apply the kinetic operator on the wavefunction, we need to Fourier Transform the wavefunction, multiply it by $\frac{\hbar^2k^2}{2m}$, and then apply the Inverse Fourier Transform. This method is computationally efficient, since the Fourier Transform can be done using the Fast Fourier Transform algorithm, which has a complexity of $\mathcal{O}(N\log N)$, where $N$ is the number of grid points (contrary to the straightforward Discrete Fourier Transform wit $\mathcal{O}(N^2)$ complexity).

A remaining challenge is the non-commutativity of the potential and kinetic operators. Because $\hat{\mathbf{V}}$ and $\hat{\mathbf{T}}$ do not commute we cannot factorize the propagator as

\begin{equation}
\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}t}\neq\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{V}}t}\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{T}}t}.
\end{equation}

If such a factorization were possible, we could apply the exponential of $\hat{\mathbf{V}}$ in the position space and the exponential of $\hat{\mathbf{T}}$ momentum space separately. To address this, we split the total propagation time $t$ into $n$ short steps $\Delta t$, so the full propagator becomes a product of short-time propagators as

\begin{equation}
\hat{\mathbf{U}}(t)=\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}t}=\mathbf{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}\Delta t}\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}\Delta t} \cdots\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}\Delta t},
\end{equation}

where $\Delta t$ is a fixed time step. Now, each individual short-time propagator acts on the wave function and causes a small evolution. We can then use the symmetric splitting approximation

\begin{equation}
\mathrm{e}^{-\frac{\mathrm{i}}{\hbar}\hat{\mathbf{H}}\Delta t}=\mathrm{e}^{-\frac{\mathrm{i}\Delta t}{2\hbar}\hat{\mathbf{V}}}\mathrm{e}^{-\frac{\mathrm{i}\Delta t}{\hbar}\hat{\mathbf{T}}}\mathrm{e}^{-\frac{\mathrm{i}\Delta t}{2\hbar}\hat{\mathbf{V}}(x)}+\mathcal{O}(\Delta t^3),
\end{equation}

which is valid for small $\Delta t$. Now we have everything we need to propagate the wavefunction in time. The formula for propagating the wavefunction from time $t$ to time $t+\Delta t$ is

\begin{equation}
\ket{\Psi\left(x,t+\Delta t\right)}=\mathrm{e}^{-\frac{\mathrm{i}\Delta t}{2\hbar}\mathbf{V}}\mathcal{F}^{-1}\left\lbrace\mathrm{e}^{-\frac{\mathrm{i}\Delta t}{\hbar}\mathbf{T}}\mathcal{F}\left\lbrace\mathrm{e}^{-\frac{\mathrm{i}\Delta t}{2\hbar}\mathbf{V}}\ket{\Psi\left(x,t\right)}\right\rbrace\right\rbrace,
\end{equation}

where $\mathbf{T}$ is the kinetic matrix in the momentum space, and $\mathbf{V}$ is the potential matrix in the position space. The matrix exponential is now simply a mathematical problem and is described [here](matrixexponential.html#matrix-exponential).

Now that we have the wavefunction an any given time, we can calculate observables, such as the density matrix, position, momentum, etc. The expectation value of position is trivial to calculate, since the wavefunction is already in the position space. The expectation value of momentum can be calculated by Fourier Transforming the wavefunction to the momentum space and multiplying by the momentum operator. The density matrix can be calculated by taking the outer product of the wavefunction with itself.

### The Adiabatic Transform

Sometimes, we would like to see the results from the dynamics in the adiabatic basis (i.e., the basis where the potential matrix $\mathbf{V}$ is diagonal). Most of the quantum chemistry software for molecular dynamics use the adiabatic basis, since the diabatic basis is not known. For that reason, most of the surface hopping algorithms also work in the adiabatic basis so it is convenient to know how to transform the wavefunction from the diabatic basis to the adiabatic basis to compare the results. To find the transformation from the diabatic basis to the adiabatic basis, we need to find a matrix $\mathbf{U}$ that diagonalizes the potential matrix $\mathbf{V}$. To do that, we solve the eigenvalue problem

\begin{equation}
\mathbf{V}\mathbf{U}=\mathbf{U}\mathbf{E},
\end{equation}

for each coordinate of the grid, where $\mathbf{E}$ is a diagonal matrix with the eigenvalues of $\mathbf{V}$ on the diagonal. The columns of $\mathbf{U}$ are the eigenvectors of $\mathbf{V}$. The matrix $\mathbf{U}$ is the transformation matrix from the diabatic basis to the adiabatic basis. To transform the wavefunction from the diabatic basis to the adiabatic basis, we multiply the wavefunction by the transformation matrix $\mathbf{U}^\dagger$ as

\begin{equation}
\ket{\Psi_{\text{adiabatic}}\left(x,t\right)}=\mathbf{U}^\dagger\ket{\Psi_{\text{diabatic}}\left(x,t\right)}.
\end{equation}

The transformation matrix $\mathbf{U}$ can be used to transform matrix representations of any operator from the diabatic basis to the adiabatic basis.

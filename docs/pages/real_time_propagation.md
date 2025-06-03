---
title: Real Time Propagation
parent: Time Evolution in Quantum Mechanics
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Real Time Propagation<!--\label{sec:real_time_propagation}-->

Real time propagation is a computational approach used to study the time evolution of quantum systems according to the Time-Dependent Schrödinger Equation. It is essential for modeling dynamic processes such as electronic excitations, molecular vibrations, and non equilibrium phenomena in quantum chemistry. By solving the wavefunction’s evolution step by step in real time, this method captures transient states and ultrafast reactions. However, due to the exponential growth of computational complexity, real time propagation is feasible only for small systems or simplified models.

## Theoretical Background

The time evolution in quantum mechanics is governed by the Time-Dependent Schrödinger Equation in the form

$$
\begin{equation}\label{eq:tdse}
\symrm{i}\hbar\frac{\symrm{d}}{\symrm{d}t}\ket{\Psi\left(t\right)}=\hat{\symbf{H}}\ket{\Psi\left(t\right)},
\end{equation}
$$

where $$i$$ is the imaginary unit, $$\hbar$$ is the reduced Planck constant, $$\ket{\Psi\left(t\right)}$$ is the time dependent wavefunction, $$\hat{\symbf{H}}$$ is the Hamiltonian operator, which encodes the total energy of the system. In the position representation, one could write $$\Psi\left(x,t\right)$$, where $$x$$ is the position coordinate for a quantum particle such as an electron or proton.

A general solution of the Time-Dependent Schrödinger Equation can be written as

$$
\begin{equation}\label{eq:tdse_propagator}
\ket{\Psi\left(t\right)}=\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}t}\ket{\Psi\left(0\right)},
\end{equation}
$$

The action of $$\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}t}$$ on the initial wavefunction $$\ket{\Psi\left(0\right)}$$ propagates the state from time $$0$$ to time $$t$$, carrying all the information about the time evolution of the system. The propagator is straightforward to find if the Hamiltonian $$\hat{\symbf{H}}$$ does not depend on time. In this case, the formal solution to the Time-Dependent Schrödinger Equation can be written as This exponential operator contains the complete description of how the system evolves with time. In practical applications, however, evaluating $$\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}t}$$ analytically is rarely straightforward for nontrivial Hamiltonians, which often necessitates numerical methods for an accurate treatment of quantum dynamics.

### Space Discretization and the Diabatic Basis

To represent the wavefunction numerically, we sample it on a spatial grid. In practical terms, this means introducing a dependence on the coordinate $$x$$, so that $$\ket{\Psi\left(t\right)}$$ becomes $$\ket{\Psi\left(x,t\right)}$$. More importantly, we expand the wavefunction $$\ket{\Psi\left(x,t\right)}$$ in a set of basis functions $$\lbrace\Phi_i\rbrace_{i=1}^N$$, where $$N$$ is the number of states included in the dynamics. In the following derivations, we will consider $$N=2$$, but the extension to multiple states (or even one state) is straightforward. The basis set will be orthonormal and satisfy the condition

$$
\begin{equation}\label{eq:diabatic_condition}
\bra{\Phi_i}\frac{\partial}{\partial x}\ket{\Phi_j}=0,\,i\neq j.
\end{equation}
$$

Because these diabatic basis states are frequently taken simply as the canonical basis in $$\mathbb{R}^N$$ when defining model potential energy surfaces, their explicit functional forms are not crucial in most model-based simulations. Nonetheless, from a physical standpoint, these states are usually determined through quantum chemical computations that ensure minimal coordinate dependence in the electronic mixing. For a two-state system, the wavefunction in the diabatic basis can be written as

$$
\begin{equation}\label{eq:diabatic_wavefunction}
\ket{\Psi\left(x,t\right)}=c_1\left(x,t\right)\ket{\Phi_1}+c_2\left(x,t\right)\ket{\Phi_2}=
\begin{pmatrix}
c_1\left(x,t\right) \\
c_2\left(x,t\right)
\end{pmatrix},
\end{equation}
$$

where $$c_i\left(x,t\right)$$ are the expansion coefficients of the wavefunction in the diabatic basis. The Hamiltonian operator $$\hat{\symbf{H}}$$ in the diabatic basis takes the form

$$
\begin{equation}\label{eq:diabatic_hamiltonian}
\hat{\symbf{H}}=\hat{\symbf{T}}+\hat{\symbf{V}}=-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}\symbf{I}+
\begin{pmatrix}
V_{11}(x) & V_{12}(x) \\
V_{21}(x) & V_{22}(x)
\end{pmatrix},
\end{equation}
$$

where $$\symbf{I}$$ is the identity matrix, $$\hat{\symbf{T}}$$ is the kinetic energy operator, and $$\symbf{V}$$ is the potential energy matrix. The diagonal elements $$V_{ii}(x)$$ are the potential energy surfaces for the diabatic state $$i$$, while the off-diagonal elements are coupling between the two states.

### Wavefunction Propagation

We aim to reduce the action of the propagator on the diabatic wavefunction \eqref{eq:diabatic_wavefunction} to a simple matrix multiplication. The main difficulty is the differential form of the kinetic operator in the equation \eqref{eq:diabatic_hamiltonian}. To handle this, we use a Fourier-based method for applying the kinetic operator $$\hat{\symbf{T}}$$ on a wavefunction $$\ket{\Psi\left(x,t\right)}$$. Taking the Fourier Transform of $$\hat{\symbf{T}}\ket{\Psi\left(x,t\right)}$$ yields

$$
\begin{align}\label{eq:kinetic_fourier_method}
\mathcal{F}\left\lbrace\hat{\symbf{T}}\ket{\Psi\left(x,t\right)}\right\rbrace&=\mathcal{F}\left\lbrace-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}\ket{\Psi\left(x,t\right)}\right\rbrace=\frac{\hbar^2 k^2}{2m}\mathcal{F}\left\lbrace\ket{\Psi\left(x,t\right)}\right\rbrace \\
\hat{\symbf{T}}\ket{\Psi\left(x,t\right)}&=\frac{\hbar^2}{2m}\mathcal{F}^{-1}\left\lbrace k^2\mathcal{F}\left\lbrace\ket{\Psi\left(x,t\right)}\right\rbrace\right\rbrace,
\end{align}
$$

where $$k$$ is the coordinate in the Fourier space. In other words, to apply the kinetic operator on the wavefunction, we need to Fourier Transform the wavefunction, multiply it by $$\frac{\hbar^2k^2}{2m}$$, and then apply the Inverse Fourier Transform. This method is computationally efficient, since the Fourier Transform can be done using the Fast Fourier Transform algorithm, which has a complexity of $$\mathcal{O}(N\log N)$$, where $$N$$ is the number of grid points (contrary to the straightforward Discrete Fourier Transform wit $$\mathcal{O}(N^2)$$ complexity).

A remaining challenge is the non-commutativity of the potential and kinetic operators. Because $$\hat{\symbf{V}}$$ and $$\hat{\symbf{T}}$$ do not commute we cannot factorize the propagator as

$$
\begin{equation}
\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}t}\neq\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{V}}t}\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{T}}t}.
\end{equation}
$$

If such a factorization were possible, we could apply the exponential of $$\hat{\symbf{V}}$$ in the position space and the exponential of $$\hat{\symbf{T}}$$ momentum space separately. To address this, we split the total propagation time $$t$$ into $$n$$ short steps $$\Delta t$$, so the full propagator becomes a product of short-time propagators as

$$
\begin{equation}
\hat{\symbf{U}}(t)=\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}t}=\symbf{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}\Delta t}\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}\Delta t} \cdots\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}\Delta t},
\end{equation}
$$

where $$\Delta t$$ is a fixed time step. Now, each individual short-time propagator acts on the wave function and causes a small evolution. We can then use the symmetric splitting approximation

$$
\begin{equation}
\symrm{e}^{-\frac{\symrm{i}}{\hbar}\hat{\symbf{H}}\Delta t}=\symrm{e}^{-\frac{\symrm{i}\Delta t}{2\hbar}\hat{\symbf{V}}}\symrm{e}^{-\frac{\symrm{i}\Delta t}{\hbar}\hat{\symbf{T}}}\symrm{e}^{-\frac{\symrm{i}\Delta t}{2\hbar}\hat{\symbf{V}}}+\mathcal{O}(\Delta t^3),
\end{equation}
$$

which is valid for small $$\Delta t$$. Now we have everything we need to propagate the wavefunction in time. The formula for propagating the wavefunction from time $$t$$ to time $$t+\Delta t$$ is

$$
\begin{equation}
\ket{\Psi\left(x,t+\Delta t\right)}=\symrm{e}^{-\frac{\symrm{i}\Delta t}{2\hbar}\symbf{V}}\mathcal{F}^{-1}\left\lbrace\symrm{e}^{-\frac{\symrm{i}\Delta t}{\hbar}\symbf{T}}\mathcal{F}\left\lbrace\symrm{e}^{-\frac{\symrm{i}\Delta t}{2\hbar}\symbf{V}}\ket{\Psi\left(x,t\right)}\right\rbrace\right\rbrace,
\end{equation}
$$

where $$\symbf{T}$$ is the kinetic matrix in the momentum space, and $$\symbf{V}$$ is the potential matrix in the position space. The matrix exponential is now simply a mathematical problem and is described [here](me.html#matrix-exponential).

Now that we have the wavefunction an any given time, we can calculate observables, such as the density matrix, position, momentum, etc. The expectation value of position is trivial to calculate, since the wavefunction is already in the position space. The expectation value of momentum can be calculated by applying Fourier Transform on the wavefunction and multiplying by the momentum operator. The density matrix can be calculated by taking the outer product of the wavefunction with itself.

### The Adiabatic Transform

Sometimes, we would like to see the results from the dynamics in the adiabatic basis (i.e., the basis where the potential matrix $$\symbf{V}$$ is diagonal). Most of the quantum chemistry software for molecular dynamics use the adiabatic basis, since the diabatic basis is not known. For that reason, most of the surface hopping algorithms also work in the adiabatic basis so it is convenient to know how to transform the wavefunction from the diabatic basis to the adiabatic basis to compare the results. To find the transformation from the diabatic basis to the adiabatic basis, we need to find a matrix $$\symbf{U}$$ that diagonalizes the potential matrix $$\symbf{V}$$. To do that, we solve the eigenvalue problem

$$
\begin{equation}
\symbf{V}\symbf{U}=\symbf{U}\symbf{E},
\end{equation}
$$

for each coordinate of the grid, where $$\symbf{E}$$ is a diagonal matrix with the eigenvalues of $$\symbf{V}$$ on the diagonal. The columns of $$\symbf{U}$$ are the eigenvectors of $$\symbf{V}$$. The matrix $$\symbf{U}$$ is the transformation matrix from the diabatic basis to the adiabatic basis. To transform the wavefunction from the diabatic basis to the adiabatic basis, we multiply the wavefunction by the transformation matrix $$\symbf{U}^\dagger$$ as

$$
\begin{equation}
\ket{\Psi_{\text{adiabatic}}\left(x,t\right)}=\symbf{U}^\dagger\ket{\Psi_{\text{diabatic}}\left(x,t\right)}.
\end{equation}
$$

The transformation matrix $$\symbf{U}$$ can be used to transform matrix representations of any operator from the diabatic basis to the adiabatic basis.

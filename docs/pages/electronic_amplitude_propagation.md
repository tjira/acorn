---
title: Quantum Amplitude Propagation
parent: Mixed Quantum-Classical Dynamics
has_children: true
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Quantum Amplitude Propagation in Mixed Schemes<!--\label{sec:quantum_amplitude_propagation_in_mixed_themes}-->

Under the Born--Oppenheimer Approximation, the total molecular wavefunction, $$\Psi(\symbf{r},\symbf{R},t)$$, is formally separable into electronic and nuclear parts because the nuclear masses greatly exceed those of the electrons. However, when two or more electronic states approach degeneracy, the Born--Oppenheimer Approximation separation fails and nonadiabatic effects become significant. Trajectory Surface Hopping is a general framework that treats nuclei and electrons separately: nuclei are propagated classically on an adiabatic Potential Energy Surface, while electrons are treated quantum mechanically.

In Trajectory Surface Hopping, the nuclei follow Newton's equations of motion on a single adiabatic Potential Energy Surface at any instant:

$$
\begin{equation}
\symbf{M}\ddot{\symbf{R}}(t)=-\nabla_{\symbf{R}}E_k(\symbf{R}(t)),
\end{equation}
$$

where $$\symbf{M}$$ is the diagonal matrix with masses for each coordinate, $$\symbf{R}$$ is nuclear geometry, and $$E_k(\symbf{R}(t))$$ is the electronic energy of the adiabatic state $$k$$ at the nuclear configuration $$\symbf{R}(t)$$. Meanwhile, the electrons evolve according to the Time-Dependent Schrödinger Equation, with the electronic wavefunction $$\Phi(\symbf{r},t;\symbf{R}(t))$$ expanded in the instantaneous adiabatic eigenvectors of the electronic Hamiltonian $$\hat{H}^{\symrm{el}}(\symbf{r};\symbf{R}(t))$$:

$$
\begin{equation}
\Phi(\symbf{r},t;\symbf{R}(t))=\sum_k c_k(t)\phi_k(r;\symbf{R}(t)),
\end{equation}
$$

where $$\phi_k(r;\symbf{R}(t))$$ satisfies the Time-Independent Schrödinger Equation

$$
\begin{equation}
\hat{H}^{\symrm{el}}(\symbf{r};\symbf{R}(t))\phi_k(r;\symbf{R}(t))=E_k(\symbf{R}(t))\phi_k(r;\symbf{R}(t)),
\end{equation}
$$

and the complex coefficients $$c_k$$ encode the instantaneous probability amplitudes for occupying each adiabatic state. For notational brevity, we will often suppress explicit dependence on the variables where no ambiguity arises.

Substituting this ansatz into the full electronic Time-Dependent Schrödinger Equation,

$$
\begin{equation}
i\hbar\frac{\partial}{\partial t}\Phi(\symbf{r},t;\symbf{R}(t))=\hat{H}^{\symrm{el}}(\symbf{r};\symbf{R}(t))\Phi(\symbf{r},t;\symbf{R}(t)),
\end{equation}
$$

and projecting onto a particular adiabatic state $$\phi_j(r;\symbf{R}(t))$$, we obtain a set of coupled equations for the coefficients $$c_j$$:

$$
\begin{equation}\label{eq:tsh_eom}
i\hbar\frac{\symrm{d}c_j}{\symrm{d}t}=\sum_k\left(H_{jk}^{\symrm{el}}-i\hbar\Braket{\phi_j\\|\frac{\partial\phi_k}{\partial t}}\right)c_k,
\end{equation}
$$

where $$\sigma_{jk}=\Braket{\phi_j\\|\frac{\partial\phi_k}{\partial t}}$$ is the Time Derivative Coupling. In the purely adiabatic basis, $$\hat{H}^{\symrm{el}}$$ is diagonal, so electronic population transfer between states is mediated entirely by $$\sigma_{jk}$$:

$$
\begin{equation}
i\hbar\frac{\partial c_j}{\partial t}=E_j c_j-i\hbar\sum_k\sigma_{jk}c_k.
\end{equation}
$$

To evaluate $$\sigma_{jk}$$ in practice, note that each adiabatic eigenfunction depends parametrically on the instantaneous nuclear coordinates $$\symbf{R}(t)$$. By the chain rule,

$$
\begin{equation}
\sigma_{jk}=\Braket{\phi_j\\|\frac{\partial\phi_k}{\partial t}}=\dot{\symbf{R}}\cdot\Braket{\phi_j\\|\frac{\partial\phi_k}{\partial\symbf{R}}}=\symbf{v}\cdot\symbf{d}_{jk},
\end{equation}
$$

where $$\symbf{v}=\dot{\symbf{R}}$$ is the nuclear velocity and $$\symbf{d}_{jk}$$ is the Nonadiabatic Coupling Vector. Although many quantum‐chemistry packages provide Nonadiabatic Coupling Vector and thus allow computation of the Time Derivative Coupling via $$\sigma_{jk}=\symbf{v}\cdot\symbf{d}_{jk}$$, Nonadiabatic Coupling Vector can become ill-defined or numerically unstable near conical intersections (e.g., in multi-reference methods). In cases where a reliable Nonadiabatic Coupling Vector is not available, one may instead employ a direct Time Derivative Coupling approximation using wavefunction overlaps.

{:.no_toc .text-delta}

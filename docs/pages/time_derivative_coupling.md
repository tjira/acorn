---
title: The Time Derivative Coupling
parent: Mixed Quantum-Classical Dynamics
layout: default
nav_order: 2
---
{% include mathjax.html %}

# The Time Derivative Coupling<!--\label{sec:time_derivative_coupling}-->

Time Derivative Coupling is a pivotal quantity in Trajectory Surface Hopping simulations because it enters the electronic equations of motion given in Eq. [here](trajectory_surface_hopping.html#mjx-eqn:eq:tsh_eom)<!--\eqref{eq:tsh_eom}-->. In *ab initio* molecular dynamics, the Time Derivative Coupling is commonly evaluated using the full Nonadiabatic Coupling Vector as $$\sigma_{jk}=\mathbf{v}\cdot\mathbf{d}_{jk}$$, where $$\mathbf{v}$$ is the nuclear velocity. When the dynamics evolves on analytic Potential Energy Surface or the Nonadiabatic Coupling Vector is for any reason not available, it is often preferable to calculate the Time Derivative Coupling directly. The most straightforward strategy applies a finite-difference formula to the adiabatic wavefunctions

$$
\begin{equation}\label{eq:tdc_fd1}
\sigma_{jk}(t)=\frac{1}{\Delta t}\left(\Braket{\phi_j(t)\\|\phi_k(t+\Delta t)}-\Braket{\phi_j(t)\\|\phi_k(t)}\right)=\frac{1}{\Delta t}\left(\Braket{\phi_j(t)\\|\phi_k(t+\Delta t)}-\delta_{jk}\right),
\end{equation}
$$

with $$\Delta t$$ a small time step. When analytic Potential Energy Surface is available, the finite-difference estimate of the Time Derivative Coupling is straightforward to implement. In *ab initio* molecular dynamics, however, it is usually more practical to recover this quantity from the non-adiabatic coupling vector, which most electronic-structure packages provide natively. The naïve first-order finite-difference formula suffers from well-known deficiencies: its truncation error scales linearly with the time step and it is highly sensitive to numerical noise. A simple remedy is to adopt the second-order approximation

$$
\begin{equation}\label{eq:tdc_fd2}
\sigma_{jk}(t)=\frac{1}{2\Delta t}\left(\Braket{\phi_j(t)\\|\phi_k(t+\Delta t)}-\Braket{\phi_j(t)\\|\phi_k(t-\Delta t)}\right).
\end{equation}
$$

Although nominally second-order accurate, the approximation still retains significant shortcomings. In particular, it still violates the fundamental antisymmetry requirement, $$\sigma_{jk}\neq-\sigma_{kj}$$, undermining the Hermiticity of the electronic Hamiltonian. The following sections tackle these deficiencies and introduce more accurate, numerically robust strategies for evaluating the Time Derivative Coupling.

## The Hammes-Schiffer Tully Scheme<!--\label{sec:hammes_schiffer_tully}-->

The Hammes-Schiffer Tully scheme remedies the antisymmetry defect of Eqs. \eqref{eq:tdc_fd1} and \eqref{eq:tdc_fd2}, making the electronic Hamiltonian Hermitian. To that end we first introduce a midpoint approximation for each adiabatic eigenstate $$\phi_j$$, obtained by linear interpolation between its values at the bracketing time slices

$$
\begin{equation}\label{eq:hst_wfn}
\phi_j\left(t+\frac{\Delta t}{2}\right)=\frac{1}{2}\left(\phi_j(t+\Delta t)+\phi_j(t)\right).
\end{equation}
$$

This midpoint state is not strictly normalised, a fully norm-preserving variant will be introduced in the next section. The time derivative at the same midpoint is obtained by combining a backward difference for $$\phi_j(t+\Delta t)$$ with a forward difference for $$\phi_j(t)$$, yielding

$$
\begin{equation}\label{eq:hst_dwfn}
\frac{\partial \phi_j}{\partial t}\left(t+\frac{\Delta t}{2}\right)=\frac{1}{\Delta t}\left(\phi_j(t+\Delta t)-\phi_j(t)\right).
\end{equation}
$$

The Time Derivative Coupling $$\sigma_{jk}=\Braket{\phi_j\\|\frac{\partial \phi_k}{\partial t}}$$ at a time $$t+\frac{\Delta t}{2}$$ is then formed as an overlap beween the expressions in Eq. \eqref{eq:hst_wfn} and Eq. \eqref{eq:hst_dwfn}:

$$
\begin{equation}\label{eq:hst_tdc}
\sigma_{jk}\left(t+\frac{\Delta t}{2}\right)=\frac{1}{2\Delta t}\left(\Braket{\phi_j(t+\Delta t)\\|\phi_k(t)}-\Braket{\phi_j(t+\Delta t)\\|\phi_k(t)}\right).
\end{equation}
$$

The Hammes-Schiffer Tully scheme in Eq. \eqref{eq:hst_tdc} is a finite-difference approximation of the Time Derivative Coupling, which is also antisymmetric, i.e., $$\sigma_{jk}=-\sigma_{kj}$$. The scheme is numerically stable and robust, but it does not preserve the norm of the adiabatic wavefunction. Any rapidly varying wavefunction will suffer from a significant normalization error, which can lead to unphysical results. To address this issue, we can apply a Norm-Preserving Interpolation to the adiabatic wavefunction, which is discussed in the next section.

## The Norm-Preserving Interpolation<!--\label{sec:norm_preserving_interpolation}-->

## The Baeck--An Scheme<!--\label{sec:baeck_an}-->

{:.no_toc .text-delta}

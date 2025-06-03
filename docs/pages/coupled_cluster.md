---
title: Coupled Cluster Theory
parent: Electronic Structure Methods
layout: default
nav_order: 4
---
{% include mathjax.html %}

# Coupled Cluster Theory

Coupled Cluster theory is a post-Hartree--Fock method used in quantum chemistry to achieve highly accurate solutions to the electronic Schrödinger equation, particularly for ground states and certain excited states. It improves upon Hartree--Fock by incorporating electron correlation effects through a systematic inclusion of excitations (singles, doubles, triples, etc.) from a reference wavefunction, usually the Hartree--Fock wavefunction. The method uses an exponential ansatz to account for these excitations, leading to a size-consistent and size-extensive approach, making it one of the most accurate methods available for small to medium-sized molecular systems.

Within Coupled Cluster theory, specific truncations are often applied to manage computational cost. The Coupled Cluster Doubles method considers only double excitations, capturing electron correlation more effectively than simpler methods like Hartree--Fock, but at a lower computational expense than higher-level methods. Coupled Cluster Singles and Doubles extends this approach by including both single and double excitations, offering greater accuracy, particularly for systems where single excitations play a significant role. Coupled Cluster Singles and Doubles is widely used due to its balance between accuracy and computational feasibility, making it a reliable choice for many chemical systems.

## Coupled Cluster Formalism

In the Coupled Cluster formalism, we write the total wavefunction in an exponential form as

$$
\begin{equation}
\ket{\Psi}=e^{\hat{\symbf{T}}}\ket{\Psi_0}
\end{equation}
$$

where $$\ket{\Psi_0}$$ is the reference wavefunction, usually the Hartree--Fock wavefunction, and $$\hat{T}$$ is the cluster operator that generates excitations from the reference wavefunction. The cluster operator is defined as

$$
\begin{equation}
\hat{\symbf{T}}=\hat{\symbf{T}}_1+\hat{\symbf{T}}_2+\hat{\symbf{T}}_3+\dots
\end{equation}
$$

where $$\hat{\symbf{T}}_1$$ generates single excitations, $$\hat{\symbf{T}}_2$$ generates double excitations, and so on. For example

$$
\begin{equation}
\hat{\symbf{T}}_1\ket{\Psi_0}=\left(\frac{1}{1!}\right)^2t_i^a\ket{\Psi_i^a}
\end{equation}
$$

where $$t_i^a$$ are the single excitation amplitudes. These amplitudes are just expansion coefficients that determine the contribution of each excitation to the total wavefunction. In the context of configuration interaction, we denoted these coefficients as $$c_i^a$$. Now that we have the total wavefunction, we want to solve the Schrödinger equation

\begin{equation}
\hat{\symbf{H}}\ket{\Psi}=E\ket{\Psi}
\end{equation}

where $$\hat{H}$$ is the molecular Hamiltonian operator, $$E$$ is the total energy of the system, and $$\ket{\Psi}$$ is the total wavefunction. In the Coupled Cluster theory, we usually rewrite the Schrödinger equation in the exponential form as

$$
\begin{equation}
e^{-\hat{\symbf{T}}}\hat{\symbf{H}}e^{\hat{\symbf{T}}}\ket{\Psi_0}=E\ket{\Psi_0}
\end{equation}
$$

because we can then express the Coupled Cluster energy as

$$
\begin{equation}
E=\braket{\Psi_0\vert e^{-\hat{\symbf{T}}}\hat{\symbf{H}}e^{\hat{\symbf{T}}}\vert\Psi_0}
\end{equation}
$$

taking advantage of the exponential form of the wavefunction. We could then proceed to express the total energy for various Coupled Cluster methods like Coupled Cluster Doubles and Coupled Cluster Singles and Doubles, but the equations would be quite lengthy. Instead, we will leave the theory here and proceed to the actual calculations. One thing to keep in mind is that the Coupled Cluster equations are nonlinear and require iterative solution methods to obtain the final amplitudes.

## Implementation of Truncated Coupled Cluster Methods

We will not go into the details here, but we will provide the final expressions for the Coupled Cluster Doubles and Coupled Cluster Singles and Doubles methods.<!--\supercite{10.1063/1.460620}--> The Coupled Cluster Doubles and Coupled Cluster Singles and Doubles methods are the most commonly used Coupled Cluster methods, and they are often used as benchmarks for other methods. All we need for the evaluation of the expressions below are the two-electron integrals in the Molecular Spinorbital basis and physicists' notation, Fock matrix in the Molecular Spinorbital basis and the orbital energy tensors obtained from the Hartree--Fock calculation. All these transformations are already explained [here](hartree_fock.html#integral-transforms-to-the-basis-of-molecular-spinorbitals)<!--in Section \ref{sec:integral_transform}-->. The expressions for the Coupled Cluster Doubles can be written as

$$
\begin{equation}
E_{\text{CCD}}=\frac{1}{4}\braket{ij\vert\vert ab}t_{ij}^{ab}
\end{equation}
$$

where the double excitation amplitudes $$t_{ij}^{ab}$$ are determined by solving the Coupled Cluster Doubles amplitude equation. The Coupled Cluster Doubles amplitude equations are given by

$$
\begin{align}
t_{ij}^{ab}=&\braket{ab\vert\vert ij}+\frac{1}{2}\braket{ab\vert\vert cd}t_{cd}^{ij}+\frac{1}{2}\braket{kl\vert\vert ij}t_{ab}^{kl}+\hat{P}_{(a/b)}\hat{P}_{(i/j)}\braket{ak\vert\vert ic}t_{cb}^{ij}-\nonumber \\
&-\frac{1}{2}\hat{P}_{(a/b)}\braket{kl\vert\vert cd}t_{ac}^{ij}t_{bd}^{kl}-\frac{1}{2}\hat{P}_{(i/j)}\braket{kl\vert\vert cd}t_{ab}^{ik}t_{cd}^{jl}+\nonumber \\
&+\frac{1}{4}\braket{kl\vert\vert cd}t_{cd}^{ij}t_{ab}^{kl}+\hat{P}_{(i/j)}\braket{kl\vert\vert cd}t_{ac}^{ik}t_{bd}^{jl}
\end{align}
$$

where $$\hat{P}_{(a/b)}$$ and $$\hat{P}_{(i/j)}$$ are permutation operators that ensure the correct antisymmetry of the amplitudes. The Coupled Cluster Singles and Doubles energy expression is given by

$$
\begin{equation}
E_{\text{CCSD}}=F_{ia}^{\symrm{MS}}t_a^i+\frac{1}{4}\braket{ij\vert\vert ab}t_{ij}^{ab}+\frac{1}{2}\braket{ij\vert\vert ab}t_{i}^{a}t_{b}^{j}
\end{equation}
$$

where the single and double excitation amplitudes $$t_a^i$$ and $$t_{ij}^{ab}$$ are determined by solving the Coupled Cluster Singles and Doubles amplitude equations. To simplify the notation a little bit, we define the the $$\mathscr{F}$$ and $$\mathscr{W}$$ intermediates as

$$
\begin{align}
\mathscr{F}_{ae}=&\left(1-\delta_{ae}\right)F_{ae}-\frac{1}{2}F_{me}t_m^a+t_m^f\braket{ma\vert\vert fe}-\frac{1}{2}\tilde{\tau}_{mn}^{af}\braket{mn\vert\vert ef} \\
\mathscr{F}_{mi}=&\left(1-\delta_{mi}\right)F_{mi}+\frac{1}{2}F_{me}t_i^e+t_n^e\braket{mn\vert\vert ie}+\frac{1}{2}\tilde{\tau}_{in}^{ef}\braket{mn\vert\vert ef} \\
\mathscr{F}_{me}=&F_{me}+t_n^f\braket{mn\vert\vert ef} \\
\mathscr{W}_{mnij}=&\braket{mn\vert\vert ij}+\hat{P}_{(i/j)}t_j^e\braket{mn\vert\vert ie}+\frac{1}{4}\tau_{ij}^{ef}\braket{mn\vert\vert ef} \\
\mathscr{W}_{abef}=&\braket{ab\vert\vert ef}-\hat{P}_{(a/b)}t_m^b\braket{am\vert\vert ef}+\frac{1}{4}\tau_{mn}^{ab}\braket{mn\vert\vert ef} \\
\mathscr{W}_{mbej}=&\braket{mb\vert\vert ej}+t_j^f\braket{mb\vert\vert ef}-t_n^b\braket{mn\vert\vert ej}-\left(\frac{1}{2}t_{jn}^{fb}+t_j^ft_n^b\right)\braket{mn\vert\vert ef}
\end{align}
$$

and two-particle excitation operators as

$$
\begin{align}
\tilde{\tau}_{ij}^{ab}=&t_{ij}^{ab}+\frac{1}{2}\left(t_i^at_j^b-t_i^bt_j^a\right) \\
\tau_{ij}^{ab}=&t_{ij}^{ab}+t_i^at_j^b-t_i^bt_j^a
\end{align}
$$

The Coupled Cluster Singles and Doubles single excitations amplitude equations are then given by

$$
\begin{align}
t_i^a=&F_{ai}^{\symrm{MS}}+t_i^e\mathscr{F}_{ae}-t_m^a\mathscr{F}_{mi}t_{im}^{ae}\mathscr{F}_{me}-t_n^f\braket{na\vert\vert if}-\nonumber-\frac{1}{2}t_{im}^{ef}\braket{ma\vert\vert ef}- \\
&-\frac{1}{2}t_{mn}^{ae}\braket{nm\vert\vert ei}
\end{align}
$$

and the Coupled Cluster Singles and Doubles double excitations amplitude equations are given by

$$
\begin{align}
t_{ij}^{ab}=&\braket{ab\vert\vert ij}+\hat{P}_{(a/b)}t_{ij}^{ae}\left(\mathscr{F}_{be}-\frac{1}{2}t_m^b\mathscr{F}_{ae}\right)-\hat{P}_{(i/j)}t_{im}^{ab}\left(\mathscr{F}_{mi}+\frac{1}{2}t_j^e\mathscr{F}_{me}\right)+\nonumber \\
&+\frac{1}{2}\tau_{mn}^{ab}\mathscr{W}_{mnij}+\frac{1}{2}\tau_{ij}^{ef}\mathscr{W}_{abef}+\hat{P}_{(i/j)}\hat{P}_{(a/b)}\left(t_{im}^{ae}\mathscr{W}_{mbej}-t_i^et_m^a\braket{mb\vert\vert ej}\right)+\nonumber \\
&+\hat{P}_{(i/j)}t_i^e\braket{ab\vert\vert ej}-\hat{P}_{(a/b)}t_m^a\braket{mb\vert\vert ij}
\end{align}
$$

The Coupled Cluster Singles and Doubles amplitude equations are, again, nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

{:.cite}
> Stanton, John F., Jürgen Gauss, John D. Watts, and Rodney J. Bartlett. 1991. “A Direct Product Decomposition Approach for Symmetry Exploitation in Many-Body Methods. I. Energy Calculations.” *The Journal of Chemical Physics* 94: 4334–45. <https://doi.org/10.1063/1.460620>.
>

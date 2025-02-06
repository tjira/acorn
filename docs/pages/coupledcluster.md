---
title: Coupled Cluster Theory
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Coupled Cluster Theory

Coupled Cluster theory is a post-Hartree--Fock method used in quantum chemistry to achieve highly accurate solutions to the electronic Schrödinger equation, particularly for ground states and certain excited states. It improves upon Hartree--Fock by incorporating electron correlation effects through a systematic inclusion of excitations (singles, doubles, triples, etc.) from a reference wavefunction, usually the Hartree--Fock wavefunction. The method uses an exponential ansatz to account for these excitations, leading to a size-consistent and size-extensive approach, making it one of the most accurate methods available for small to medium-sized molecular systems.

Within Coupled Cluster theory, specific truncations are often applied to manage computational cost. The Coupled Cluster Doubles method considers only double excitations, capturing electron correlation more effectively than simpler methods like Hartree--Fock, but at a lower computational expense than higher-level methods. Coupled Cluster Singles and Doubles extends this approach by including both single and double excitations, offering greater accuracy, particularly for systems where single excitations play a significant role. Coupled Cluster Singles and Doubles is widely used due to its balance between accuracy and computational feasibility, making it a reliable choice for many chemical systems.
## Coupled Cluster Formalism

In the Coupled Cluster formalism, we write the total wavefunction in an exponential form as

\begin{equation}
\ket{\Psi}=e^{\hat{\mathbf{T}}}\ket{\Psi\_0}
\end{equation}

where $\ket{\Psi\_0}$ is the reference wavefunction, usually the Hartree--Fock wavefunction, and $\hat{T}$ is the cluster operator that generates excitations from the reference wavefunction. The cluster operator is defined as

\begin{equation}
\hat{\mathbf{T}}=\hat{\mathbf{T}}\_1+\hat{\mathbf{T}}\_2+\hat{\mathbf{T}}\_3+\dots
\end{equation}

where $\hat{\mathbf{T}}\_1$ generates single excitations, $\hat{\mathbf{T}}\_2$ generates double excitations, and so on. For example

\begin{equation}
\hat{\mathbf{T}}\_1\ket{\Psi\_0}=\left(\frac{1}{1!}\right)^2t\_i^a\ket{\Psi\_i^a}
\end{equation}

where $t\_i^a$ are the single excitation amplitudes. These amplitudes are just expansion coefficients that determine the contribution of each excitation to the total wavefunction. In the context of configuration interaction, we denoted these coefficients as $c\_i^a$. Now that we have the total wavefunction, we want to solve the Schrödinger equation

\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi}
\end{equation}

where $\hat{H}$ is the molecular Hamiltonian operator, $E$ is the total energy of the system, and $\ket{\Psi}$ is the total wavefunction. In the Coupled Cluster theory, we usually rewrite the Schrödinger equation in the exponential form as

\begin{equation}
e^{-\hat{\mathbf{T}}}\hat{\mathbf{H}}e^{\hat{\mathbf{T}}}\ket{\Psi\_0}=E\ket{\Psi\_0}
\end{equation}

because we can then express the Coupled Cluster energy as

\begin{equation}
E=\braket{\Psi\_0|e^{-\hat{\mathbf{T}}}\hat{\mathbf{H}}e^{\hat{\mathbf{T}}}|\Psi\_0}
\end{equation}

taking advantage of the exponential form of the wavefunction. We could then proceed to express the total energy for various Coupled Cluster methods like Coupled Cluster Doubles and Coupled Cluster Singles and Doubles, but the equations would be quite lengthy. Instead, we will leave the theory here and proceed to the actual calculations. One thing to keep in mind is that the Coupled Cluster equations are nonlinear and require iterative solution methods to obtain the final amplitudes.

## Implementation of Truncated Coupled Cluster Methods

We will not go into the details here, but we will provide the final expressions for the Coupled Cluster Doubles and Coupled Cluster Singles and Doubles methods.<!--\supercite{10.1063/1.460620}--> The Coupled Cluster Doubles and Coupled Cluster Singles and Doubles methods are the most commonly used Coupled Cluster methods, and they are often used as benchmarks for other methods. All we need for the evaluation of the expressions below are the two-electron integrals in the Molecular Spinorbital basis and physicists' notation, Fock matrix in the Molecular Spinorbital basis and the orbital energy tensors obtained from the Hartree--Fock calculation. All these transformations are already explained [here](hartreefockmethod.html#integral-transforms-to-the-basis-of-molecular-spinorbitals). The expressions for the Coupled Cluster Doubles can be written as

\begin{equation}
E\_{\text{CCD}}=\frac{1}{4}\braket{ij||ab}t\_{ij}^{ab}
\end{equation}

where the double excitation amplitudes $t\_{ij}^{ab}$ are determined by solving the Coupled Cluster Doubles amplitude equation. The Coupled Cluster Doubles amplitude equations are given by

\begin{align}
t\_{ij}^{ab}=&\braket{ab||ij}+\frac{1}{2}\braket{ab||cd}t\_{cd}^{ij}+\frac{1}{2}\braket{kl||ij}t\_{ab}^{kl}+\hat{P}\_{(a/b)}\hat{P}\_{(i/j)}\braket{ak||ic}t\_{cb}^{ij}-\nonumber \\\\\
&-\frac{1}{2}\hat{P}\_{(a/b)}\braket{kl||cd}t\_{ac}^{ij}t\_{bd}^{kl}-\frac{1}{2}\hat{P}\_{(i/j)}\braket{kl||cd}t\_{ab}^{ik}t\_{cd}^{jl}+\nonumber \\\\\
&+\frac{1}{4}\braket{kl||cd}t\_{cd}^{ij}t\_{ab}^{kl}+\hat{P}\_{(i/j)}\braket{kl||cd}t\_{ac}^{ik}t\_{bd}^{jl}
\end{align}

where $\hat{P}\_{(a/b)}$ and $\hat{P}\_{(i/j)}$ are permutation operators that ensure the correct antisymmetry of the amplitudes. The Coupled Cluster Singles and Doubles energy expression is given by

\begin{equation}
E\_{\text{CCSD}}=F\_{ia}^{\mathrm{MS}}t\_a^i+\frac{1}{4}\braket{ij||ab}t\_{ij}^{ab}+\frac{1}{2}\braket{ij||ab}t\_{i}^{a}t\_{b}^{j}
\end{equation}

where the single and double excitation amplitudes $t\_a^i$ and $t\_{ij}^{ab}$ are determined by solving the Coupled Cluster Singles and Doubles amplitude equations. To simplify the notation a little bit, we define the the $\mathscr{F}$ and $\mathscr{W}$ intermediates as

\begin{align}
\mathscr{F}\_{ae}=&\left(1-\delta\_{ae}\right)F\_{ae}-\frac{1}{2}F\_{me}t\_m^a+t\_m^f\braket{ma||fe}-\frac{1}{2}\tilde{\tau}\_{mn}^{af}\braket{mn||ef} \\\\\
\mathscr{F}\_{mi}=&\left(1-\delta\_{mi}\right)F\_{mi}+\frac{1}{2}F\_{me}t\_i^e+t\_n^e\braket{mn||ie}+\frac{1}{2}\tilde{\tau}\_{in}^{ef}\braket{mn||ef} \\\\\
\mathscr{F}\_{me}=&F\_{me}+t\_n^f\braket{mn||ef} \\\\\
\mathscr{W}\_{mnij}=&\braket{mn||ij}+\hat{P}\_{(i/j)}t\_j^e\braket{mn||ie}+\frac{1}{4}\tau\_{ij}^{ef}\braket{mn||ef} \\\\\
\mathscr{W}\_{abef}=&\braket{ab||ef}-\hat{P}\_{(a/b)}t\_m^b\braket{am||ef}+\frac{1}{4}\tau\_{mn}^{ab}\braket{mn||ef} \\\\\
\mathscr{W}\_{mbej}=&\braket{mb||ej}+t\_j^f\braket{mb||ef}-t\_n^b\braket{mn||ej}-\left(\frac{1}{2}t\_{jn}^{fb}+t\_j^ft\_n^b\right)\braket{mn||ef}
\end{align}

and two-particle excitation operators as

\begin{align}
\tilde{\tau}\_{ij}^{ab}=&t\_{ij}^{ab}+\frac{1}{2}\left(t\_i^at\_j^b-t\_i^bt\_j^a\right) \\\\\
\tau\_{ij}^{ab}=&t\_{ij}^{ab}+t\_i^at\_j^b-t\_i^bt\_j^a
\end{align}

The Coupled Cluster Singles and Doubles single excitations amplitude equations are then given by

\begin{align}
t\_i^a=&F\_{ai}^{\mathrm{MS}}+t\_i^e\mathscr{F}\_{ae}-t\_m^a\mathscr{F}\_{mi}t\_{im}^{ae}\mathscr{F}\_{me}-t\_n^f\braket{na||if}-\nonumber-\frac{1}{2}t\_{im}^{ef}\braket{ma||ef}- \\\\\
&-\frac{1}{2}t\_{mn}^{ae}\braket{nm||ei}
\end{align}

and the Coupled Cluster Singles and Doubles double excitations amplitude equations are given by

\begin{align}
t\_{ij}^{ab}=&\braket{ab||ij}+\hat{P}\_{(a/b)}t\_{ij}^{ae}\left(\mathscr{F}\_{be}-\frac{1}{2}t\_m^b\mathscr{F}\_{ae}\right)-\hat{P}\_{(i/j)}t\_{im}^{ab}\left(\mathscr{F}\_{mi}+\frac{1}{2}t\_j^e\mathscr{F}\_{me}\right)+\nonumber \\\\\
&+\frac{1}{2}\tau\_{mn}^{ab}\mathscr{W}\_{mnij}+\frac{1}{2}\tau\_{ij}^{ef}\mathscr{W}\_{abef}+\hat{P}\_{(i/j)}\hat{P}\_{(a/b)}\left(t_{im}^{ae}\mathscr{W}\_{mbej}-t\_i^et\_m^a\braket{mb||ej}\right)+\nonumber \\\\\
&+\hat{P}\_{(i/j)}t\_i^e\braket{ab||ej}-\hat{P}\_{(a/b)}t\_m^a\braket{mb||ij}
\end{align}

The Coupled Cluster Singles and Doubles amplitude equations are, again, nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

## Coupled Cluster Singles and Doubles Code Exercise

After completing the Hartree--Fock implementation [here](hartreefockmethod.html#hartreefock-method-and-integral-transform-coding-exercise), you can proceed with coding the Coupled Cluster Singles and Doubles exercise, which builds on the Hartree--Fock results. The exercise is provided in the Listing <!--\ref{code:cc_exercise}--> below.

<!--{id=code:cc_exercise caption="Coupled Cluster Singles and Doubles exercise code. The energy and amplitudes are initialized with default values. The student is expected to fill the loop for the calculation of the excitation amplitudes and ground state energy. After the self-consistency is achieved the result is automatically printed."}-->
```python
"""
We also have everything we need for the CC calculations. In this exercise, we will calculate the CCSD energy. Since the calculation will be iterative, I define here the CCSD energy as zero, the "E_CCSD_P" variable will be used to monitor convergence.
"""
E_CCSD, E_CCSD_P = 0, 1

"""
The first step of the calculation is to define the "t1" and "t2" amplitudes. These arrays can be initialized as zero arrays with the appropriate dimensions. I will leave this task to you.
"""
t1, t2 = np.array([]), np.array([])

"""
Now for the more complicated part. The CCSD calculation is iterative, and the convergence criterion is set by the "thresh" variable. The while loop should be filled with the appropriate calculations. The calculation of the "t1" and "t2" amplitudes is the most challenging part of the CCSD calculation. After convergence, the "E_CCSD" variable should store the final CCSD energy.
"""
while abs(E_CCSD - E_CCSD_P) > thresh:
    break

# print the CCSD energy
print("CCSD ENERGY: {:.8f}".format(E_HF + E_CCSD + VNN))
```

Solution to this exercise can be found [here](codesolutions.html#coupled-cluster-singles-and-doubles).

{:.cite}
> Stanton, John F., Jürgen Gauss, John D. Watts, and Rodney J. Bartlett. 1991. “A Direct Product Decomposition Approach for Symmetry Exploitation in Many-Body Methods. I. Energy Calculations.” *The Journal of Chemical Physics* 94: 4334–45. <https://doi.org/10.1063/1.460620>.
>

---
title: Coupled Clusters Theory
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Coupled Cluster Theory

Coupled-cluster theory is a post-Hartree-Fock method used in quantum chemistry to achieve highly accurate solutions to the electronic Schrödinger equation, particularly for ground states and certain excited states. It improves upon Hartree-Fock by incorporating electron correlation effects through a systematic inclusion of excitations (singles, doubles, triples, etc.) from a reference wavefunction, usually the Hartree-Fock wavefunction. The method uses an exponential ansatz to account for these excitations, leading to a size-consistent and size-extensive approach, making it one of the most accurate methods available for small to medium-sized molecular systems.

## Base Coupled Clusters Methods

Within coupled-cluster theory, specific truncations are often applied to manage computational cost. The Coupled-Cluster with Doubles (CCD) method considers only double excitations, capturing electron correlation more effectively than simpler methods like Hartree-Fock but at a lower computational expense than higher-level methods. Coupled-Cluster with Singles and Doubles (CCSD) extends this approach by including both single and double excitations, offering greater accuracy, particularly for systems where single excitations play a significant role. CCSD is widely used due to its balance between accuracy and computational feasibility, making it a reliable choice for many chemical systems.

### The Transformation of Integrals to Molecular Spinorbital Basis

To begin, we need to convert the Fock matrix $\mathbf{F}$ and the Coulomb integrals $\mathbf{J}$ into the basis of molecular spinorbitals (MS). The proces of transforming the Coulomb integrals is already described in the Møller–Plesset perturbation theory page [here](mollerplessetperturbationtheory.html#the-transform-of-integrals-to-molecular-spinorbital-basis). To transform the Fock matrix, we use the formula

\begin{equation}
F_{pq}^{MS}=C_{\mu p}^{MS}(\mathbf{I}\_{2}\otimes_K\mathbf{F})\_{\mu\nu}C_{\nu q}^{MS}
\end{equation}

where $\mathbf{C}^{MS}$ is the coefficient matrix in the MS basis, defined in the MPPT page. Remember that the sums over $i$, $j$, $k$, and $l$ are over occupied spinorbitals, while the sums over $a$, $b$, $c$, and $d$ are over virtual spinorbitals.

### Coupled Clusters Doubles (CCD)

The CCD energy expression is given by

\begin{equation}
E_{\text{CCD}}=\frac{1}{4}\braket{ij||ab}t_{ab}^{ij}
\end{equation}

where the double excitation amplitudes $t_{ab}^{ij}$ are determined by solving the CCD amplitude equation. The CCD amplitude equations are given by

\begin{align}
t_{ab}^{ij}=&\braket{ab||ij}+\frac{1}{2}\braket{ab||cd}t_{cd}^{ij}+\frac{1}{2}\braket{kl||ij}t_{ab}^{kl}+\hat{P}\_{(a/b)}\hat{P}\_{(i/j)}\braket{ak||ic}t_{cb}^{ij}-\frac{1}{2}\hat{P}\_{(a/b)}\braket{kl||cd}t_{ac}^{ij}t_{bd}^{kl} - \\\\\
&-\frac{1}{2}\hat{P}\_{(i/j)}\braket{kl||cd}t_{ab}^{ik}t_{cd}^{jl}+\frac{1}{4}\braket{kl||cd}t_{cd}^{ij}t_{ab}^{kl}+\hat{P}\_{(i/j)}\braket{kl||cd}t_{ac}^{ik}t_{bd}^{jl}
\end{align}

where $\hat{P}\_{(a/b)}$ and $\hat{P}\_{(i/j)}$ are permutation operators that ensure the correct antisymmetry of the amplitudes. As you can see, the CCD amplitude equations are nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

### Coupled Clusters Singles and Doubles (CCSD)

The CCSD energy expression is given by

\begin{equation}
E_{\text{CCSD}}=\mathbf{F}\_{ia}^{MS}t_a^i+\frac{1}{4}\braket{ij||ab}t_{ab}^{ij}+\frac{1}{2}\braket{ij||ab}t_{a}^{i}t_{b}^{j}
\end{equation}

where the single and double excitation amplitudes $t_a^i$ and $t_{ab}^{ij}$ are determined by solving the CCSD amplitude equations. To simplify the notation a little bit, we define the the $\mathscr{F}$ and $\mathscr{W}$ intermediates as


\begin{align}
\mathscr{F}\_{ae}=&\left(1-\delta\_{ae}\right)F\_{ae}-\frac{1}{2}\sum\_mF_{me}t_a^m+\sum\_{mf}t\_f^m\braket{ma||fe}-\frac{1}{2}\sum\_{mnf}\tilde{\tau}\_{af}^{mn}\braket{mn||ef} \\\\\
\mathscr{F}\_{mi}=&\left(1-\delta\_{mi}\right)F\_{mi}+\frac{1}{2}\sum\_eF_{me}t_e^i+\sum\_{en}t\_e^n\braket{mn||ie}+\frac{1}{2}\sum\_{nef}\tilde{\tau}\_{ef}^{in}\braket{mn||ef} \\\\\
\mathscr{F}\_{me}=&F_{me}+\sum\_{nf}t\_f^n\braket{mn||ef} \\\\\
\mathscr{W}\_{mnij}=& \\\\\
\mathscr{W}\_{abef}=& \\\\\
\mathscr{W}\_{mbej}=&
\end{align}

and two-particle excitation operators as

\begin{align}
\tilde{\tau}\_{ab}^{ij}=&t_{ab}^{ij}+\frac{1}{2}\left(t_a^it_b^j-t_b^it_a^j\right) \\\\\
\tau\_{ab}^{ij}=&t_{ab}^{ij}+t_a^it_b^j-t_b^it_a^j
\end{align}

The CCSD single amplitude equations are then given by

\begin{equation}
t_a^i=
\end{equation}

The CCSD double amplitude equations are given by

\begin{equation}
t_{ab}^{ij}=
\end{equation}

The CCSD amplitude equations are, again, nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

{:.cite}
> John F. Stanton, Jurgen Gauss, John D. Watts, and Rodney J. Bartlett 1991. Journal of Chemical Physics, 94, p.4334-4345.

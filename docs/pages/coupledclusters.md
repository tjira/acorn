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
E_{\text{CCD}}=\frac{1}{4}\braket{ij||ab}t_{ij}^{ab}
\end{equation}

where the double excitation amplitudes $t_{ij}^{ab}$ are determined by solving the CCD amplitude equation. The CCD amplitude equations are given by

\begin{align}
t_{ij}^{ab}=&\braket{ab||ij}+\frac{1}{2}\braket{ab||cd}t_{cd}^{ij}+\frac{1}{2}\braket{kl||ij}t_{ab}^{kl}+\hat{P}\_{(a/b)}\hat{P}\_{(i/j)}\braket{ak||ic}t_{cb}^{ij}-\frac{1}{2}\hat{P}\_{(a/b)}\braket{kl||cd}t_{ac}^{ij}t_{bd}^{kl}- \\\\\
&-\frac{1}{2}\hat{P}\_{(i/j)}\braket{kl||cd}t_{ab}^{ik}t_{cd}^{jl}+\frac{1}{4}\braket{kl||cd}t_{cd}^{ij}t_{ab}^{kl}+\hat{P}\_{(i/j)}\braket{kl||cd}t_{ac}^{ik}t_{bd}^{jl}
\end{align}

where $\hat{P}\_{(a/b)}$ and $\hat{P}\_{(i/j)}$ are permutation operators that ensure the correct antisymmetry of the amplitudes. As you can see, the CCD amplitude equations are nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

### Coupled Clusters Singles and Doubles (CCSD)

The CCSD energy expression is given by

\begin{equation}
E_{\text{CCSD}}=\mathbf{F}\_{ia}^{MS}t_a^i+\frac{1}{4}\braket{ij||ab}t_{ij}^{ab}+\frac{1}{2}\braket{ij||ab}t_{i}^{a}t_{b}^{j}
\end{equation}

where the single and double excitation amplitudes $t_a^i$ and $t_{ij}^{ab}$ are determined by solving the CCSD amplitude equations. To simplify the notation a little bit, we define the the $\mathscr{F}$ and $\mathscr{W}$ intermediates as


\begin{align}
\mathscr{F}\_{ae}=&\left(1-\delta\_{ae}\right)F\_{ae}-\frac{1}{2}\sum\_mF_{me}t_m^a+\sum\_{mf}t\_m^f\braket{ma||fe}-\frac{1}{2}\sum\_{mnf}\tilde{\tau}\_{mn}^{af}\braket{mn||ef} \\\\\
\mathscr{F}\_{mi}=&\left(1-\delta\_{mi}\right)F\_{mi}+\frac{1}{2}\sum\_eF_{me}t_i^e+\sum\_{en}t\_n^e\braket{mn||ie}+\frac{1}{2}\sum\_{nef}\tilde{\tau}\_{in}^{ef}\braket{mn||ef} \\\\\
\mathscr{F}\_{me}=&F_{me}+\sum\_{nf}t\_n^f\braket{mn||ef} \\\\\
\mathscr{W}\_{mnij}=&\braket{mn||ij}+\hat{P}\_{(i/j)}\sum_et_j^e\braket{mn||ie}+\frac{1}{4}\sum_{ef}\tau_{ij}^{ef}\braket{mn||ef} \\\\\
\mathscr{W}\_{abef}=&\braket{ab||ef}-\hat{P}\_{(a/b)}\sum_et_m^b\braket{am||ef}+\frac{1}{4}\sum_{mn}\tau_{mn}^{ab}\braket{mn||ef} \\\\\
\mathscr{W}\_{mbej}=&\braket{mb||ej}+\sum_ft_j^f\braket{mb||ef}-\sum_nt_n^b\braket{mn||ej}-\sum_{nf}\left(\frac{1}{2}t_{jn}^{fb}+t_j^ft_n^b\right)\braket{mn||ef}
\end{align}

and two-particle excitation operators as

\begin{align}
\tilde{\tau}\_{ij}^{ab}=&t_{ij}^{ab}+\frac{1}{2}\left(t_i^at_j^b-t_i^bt_j^a\right) \\\\\
\tau\_{ij}^{ab}=&t_{ij}^{ab}+t_i^at_j^b-t_i^bt_j^a
\end{align}

The CCSD single amplitude equations are then given by

\begin{align}
t_i^a=&F_{ai}^{MS}+\sum_et_i^e\mathscr{F}\_{ae}-\sum_mt_m^a\mathscr{F}\_{mi}\sum_{me}t_{im}^{ae}\mathscr{F}\_{me}-\sum_{nf}t_n^f\braket{na||if}- \\\\\
&-\frac{1}{2}\sum_{mef}t_{im}^{ef}\braket{ma||ef}-\frac{1}{2}\sum_{men}t_{mn}^{ae}\braket{nm||ei}
\end{align}

The CCSD double amplitude equations are given by

\begin{align}
t_{ij}^{ab}=&\braket{ab||ij}+\hat{P}\_{(a/b)}\sum_et_{ij}^{ae}\left(\mathscr{F}\_{be}-\frac{1}{2}\sum_mt_m^b\mathscr{F}\_{ae}\right)-\hat{P}\_{(i/j)}\sum_mt_{im}^{ab}\left(\mathscr{F}\_{mi}+\frac{1}{2}\sum_et_j^e\mathscr{F}\_{me}\right)+ \\\\\
&+\frac{1}{2}\sum_{mn}\tau_{mn}^{ab}\mathscr{W}\_{mnij}+\frac{1}{2}\sum_{ef}\tau_{ij}^{ef}\mathscr{W}\_{abef}+\hat{P}\_{(i/j)}\hat{P}\_{(a/b)}\sum_{me}\left(t_{im}^{ae}\mathscr{W}\_{mbej}-t_i^et_m^a\braket{mb||ej}\right)+ \\\\\
&+\hat{P}\_{(i/j)}\sum_et_i^e\braket{ab||ej}-\hat{P}\_{(a/b)}\sum_mt_m^a\braket{mb||ij}
\end{align}

The CCSD amplitude equations are, again, nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

{:.cite}
> John F. Stanton, Jurgen Gauss, John D. Watts, and Rodney J. Bartlett 1991. Journal of Chemical Physics, 94, p.4334-4345.

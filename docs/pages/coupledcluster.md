---
title: Coupled Cluster Theory
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Coupled Cluster Theory

Coupled Cluster theory is a post-Hartree--Fock method used in quantum chemistry to achieve highly accurate solutions to the electronic Schrödinger equation, particularly for ground states and certain excited states. It improves upon Hartree--Fock by incorporating electron correlation effects through a systematic inclusion of excitations (singles, doubles, triples, etc.) from a reference wavefunction, usually the Hartree--Fock wavefunction. The method uses an exponential ansatz to account for these excitations, leading to a size-consistent and size-extensive approach, making it one of the most accurate methods available for small to medium-sized molecular systems.

Within Coupled Cluster theory, specific truncations are often applied to manage computational cost. The Coupled Cluster Doubles method considers only double excitations, capturing electron correlation more effectively than simpler methods like Hartree--Fock but at a lower computational expense than higher-level methods. Coupled Cluster Singles and Doubles extends this approach by including both single and double excitations, offering greater accuracy, particularly for systems where single excitations play a significant role. Coupled Cluster Singles and Doubles is widely used due to its balance between accuracy and computational feasibility, making it a reliable choice for many chemical systems. In the below equations, we will again use the convention that the indices $i$, $j$, $k$, and $l$ run over occupied spinorbitals, while the indices $a$, $b$, $c$, and $d$ run over virtual spinorbitals.

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
\hat{\mathbf{T}}\_1\ket{\Psi\_0}=\left(\frac{1}{1!}\right)^2t\_i^a\ket{\Psi\_i^a},
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
E=\braket{\Psi\_0|e^{-\hat{\mathbf{T}}}\hat{\mathbf{H}}e^{\hat{\mathbf{T}}}|\Psi\_0},
\end{equation}

taking advantage of the exponential form of the wavefunction. We could then proceed to express the total energy for various Coupled Cluster methods like Coupled Cluster Doubles and Coupled Cluster Singles and Doubles, but the equations would be quite lengthy. Instead, we will leave the theory here and proceed to the actual calculations. One thing to keep in mind is that the CC equations are nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

## Implementation of CCD and CCSD

The derivation of the equations that are actually used to perform the calculations is quite lengthy and involves a lot of algebra. We will not go into the details here, but we will provide the final expressions for the Coupled Cluster Doubles and Coupled Cluster Singles and Doubles methods.<!--\cite{10.1063/1.460620}--> The Coupled Cluster Doubles and Coupled Cluster Singles and Doubles methods are the most commonly used Coupled Cluster methods, and they are often used as benchmarks for other methods. All we need for the evaluation of the expressions below are the Coulomb integrals in the MS basis and physicists' notation, Fock matrix in the MS basis and the orbital energies obtained from the Hartree--Fock calculation. All these transformations are already explained [here](hartreefockmethod.html#the-integral-transforms) in the Hartree--Fock section. The expressions for the Coupled Cluster Doubles can be written as

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
E\_{\text{CCSD}}=F\_{ia}^{MS}t\_a^i+\frac{1}{4}\braket{ij||ab}t\_{ij}^{ab}+\frac{1}{2}\braket{ij||ab}t\_{i}^{a}t\_{b}^{j}
\end{equation}

where the single and double excitation amplitudes $t\_a^i$ and $t\_{ij}^{ab}$ are determined by solving the Coupled Cluster Singles and Doubles amplitude equations. To simplify the notation a little bit, we define the the $\mathscr{F}$ and $\mathscr{W}$ intermediates as

\begin{align}
\mathscr{F}\_{ae}=&\left(1-\delta\_{ae}\right)F\_{ae}-\frac{1}{2}\sum\_mF\_{me}t\_m^a+\sum\_{mf}t\_m^f\braket{ma||fe}-\frac{1}{2}\sum\_{mnf}\tilde{\tau}\_{mn}^{af}\braket{mn||ef} \\\\\
\mathscr{F}\_{mi}=&\left(1-\delta\_{mi}\right)F\_{mi}+\frac{1}{2}\sum\_eF\_{me}t\_i^e+\sum\_{en}t\_n^e\braket{mn||ie}+\frac{1}{2}\sum\_{nef}\tilde{\tau}\_{in}^{ef}\braket{mn||ef} \\\\\
\mathscr{F}\_{me}=&F\_{me}+\sum\_{nf}t\_n^f\braket{mn||ef} \\\\\
\mathscr{W}\_{mnij}=&\braket{mn||ij}+\hat{P}\_{(i/j)}\sum\_et\_j^e\braket{mn||ie}+\frac{1}{4}\sum\_{ef}\tau\_{ij}^{ef}\braket{mn||ef} \\\\\
\mathscr{W}\_{abef}=&\braket{ab||ef}-\hat{P}\_{(a/b)}\sum\_et\_m^b\braket{am||ef}+\frac{1}{4}\sum\_{mn}\tau\_{mn}^{ab}\braket{mn||ef} \\\\\
\mathscr{W}\_{mbej}=&\braket{mb||ej}+\sum\_ft\_j^f\braket{mb||ef}-\sum\_nt\_n^b\braket{mn||ej}- \\\\\
&-\sum\_{nf}\left(\frac{1}{2}t\_{jn}^{fb}+t\_j^ft\_n^b\right)\braket{mn||ef}
\end{align}

and two-particle excitation operators as

\begin{align}
\tilde{\tau}\_{ij}^{ab}=&t\_{ij}^{ab}+\frac{1}{2}\left(t\_i^at\_j^b-t\_i^bt\_j^a\right) \\\\\
\tau\_{ij}^{ab}=&t\_{ij}^{ab}+t\_i^at\_j^b-t\_i^bt\_j^a
\end{align}

The Coupled Cluster Singles and Doubles single excitations amplitude equations are then given by

\begin{align}
t\_i^a=&F\_{ai}^{MS}+\sum\_et\_i^e\mathscr{F}\_{ae}-\sum\_mt\_m^a\mathscr{F}\_{mi}\sum\_{me}t\_{im}^{ae}\mathscr{F}\_{me}-\sum\_{nf}t\_n^f\braket{na||if}-\nonumber \\\\\
&-\frac{1}{2}\sum\_{mef}t\_{im}^{ef}\braket{ma||ef}-\frac{1}{2}\sum\_{men}t\_{mn}^{ae}\braket{nm||ei}
\end{align}

and the Coupled Cluster Singles and Doubles double excitations amplitude equations are given by

\begin{align}
t\_{ij}^{ab}=&\braket{ab||ij}+\hat{P}\_{(a/b)}\sum\_et\_{ij}^{ae}\left(\mathscr{F}\_{be}-\frac{1}{2}\sum\_mt\_m^b\mathscr{F}\_{ae}\right)-\nonumber \\\\\
&-\hat{P}\_{(i/j)}\sum\_mt\_{im}^{ab}\left(\mathscr{F}\_{mi}+\frac{1}{2}\sum\_et\_j^e\mathscr{F}\_{me}\right)+\frac{1}{2}\sum\_{mn}\tau\_{mn}^{ab}\mathscr{W}\_{mnij}+\nonumber \\\\\
&+\frac{1}{2}\sum\_{ef}\tau\_{ij}^{ef}\mathscr{W}\_{abef}+\hat{P}\_{(i/j)}\hat{P}\_{(a/b)}\sum_{me}\left(t_{im}^{ae}\mathscr{W}\_{mbej}-t\_i^et\_m^a\braket{mb||ej}\right)+\nonumber \\\\\
&+\hat{P}\_{(i/j)}\sum\_et\_i^e\braket{ab||ej}-\hat{P}\_{(a/b)}\sum\_mt\_m^a\braket{mb||ij}
\end{align}

The Coupled Cluster Singles and Doubles amplitude equations are, again, nonlinear and require iterative solution methods to obtain the final amplitudes. The initial guess for the amplitudes is often set to zero, and the equations are solved iteratively until convergence is achieved.

## Code Examples

### Exercise

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

### Solution

```python
# energy containers for all the CC methods
E_CCD, E_CCD_P, E_CCSD, E_CCSD_P = 0, 1, 0, 1

# initialize the first guess for the t-amplitudes as zeros
t1, t2 = np.zeros((2 * nvirt, 2 * nocc)), np.zeros((2 * nvirt, 2 * nvirt, 2 * nocc, 2 * nocc))

# CCD loop
if args.ccd:
    while abs(E_CCD - E_CCD_P) > args.threshold:

        # collect all the distinct LCCD terms
        lccd1 = 0.5 * np.einsum("abcd,cdij->abij", Jmsa[v, v, v, v], t2, optimize=True)
        lccd2 = 0.5 * np.einsum("klij,abkl->abij", Jmsa[o, o, o, o], t2, optimize=True)
        lccd3 =       np.einsum("akic,bcjk->abij", Jmsa[v, o, o, v], t2, optimize=True)

        # apply the permuation operator and add it to the corresponding LCCD terms
        lccd3 = lccd3 + lccd3.transpose(1, 0, 3, 2) - lccd3.transpose(1, 0, 2, 3) - lccd3.transpose(0, 1, 3, 2)

        # collect all the distinct first CCD terms
        ccd1 = -0.50 * np.einsum("klcd,acij,bdkl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)
        ccd2 = -0.50 * np.einsum("klcd,abik,cdjl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)
        ccd3 =  0.25 * np.einsum("klcd,cdij,abkl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)
        ccd4 =         np.einsum("klcd,acik,bdjl->abij", Jmsa[o, o, v, v], t2, t2, optimize=True)

        # apply the permuation operator and add it to the corresponding CCD terms
        ccd1, ccd2, ccd4 = ccd1 - ccd1.transpose(1, 0, 2, 3), ccd2 - ccd2.transpose(0, 1, 3, 2), ccd4 - ccd4.transpose(0, 1, 3, 2)

        # update the t-amplitudes
        t2 = Emsd * (Jmsa[v, v, o, o] + lccd1 + lccd2 + lccd3 + ccd1 + ccd2 + ccd3 + ccd4)

        # evaluate the energy
        E_CCD_P, E_CCD = E_CCD, 0.25 * np.einsum("ijab,abij", Jmsa[o, o, v, v], t2, optimize=True)

    # print the CCD energy
    print("    CCD ENERGY: {:.8f}".format(E_HF + E_CCD + VNN))

# CCSD loop
if args.ccsd:
    while abs(E_CCSD - E_CCSD_P) > args.threshold:

        # calculate the effective two-particle excitation operators
        ttau = t2 + 0.5 * np.einsum("ai,bj->abij", t1, t1, optimize=True) - 0.5 * np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(2, 3)
        tau  = t2 +       np.einsum("ai,bj->abij", t1, t1, optimize=True) -       np.einsum("ai,bj->abij", t1, t1, optimize=True).swapaxes(2, 3)

        # calculate the 2D two-particle intermediates
        Fae = (1 - np.eye(2 * nvirt)) * Fms[v, v] - 0.5 * np.einsum("me,am->ae",     Fms[o, v],        t1,   optimize=True)                                                   +       np.einsum("mafe,fm->ae",   Jmsa[o, v, v, v], t1,   optimize=True)                                                   - 0.5 * np.einsum("mnef,afmn->ae", Jmsa[o, o, v, v], ttau, optimize=True)
        Fmi = (1 - np.eye(2 * nocc )) * Fms[o, o] + 0.5 * np.einsum("me,ei->mi",     Fms[o, v],        t1,   optimize=True)                                                   +       np.einsum("mnie,en->mi",   Jmsa[o, o, o, v], t1,   optimize=True)                                                   + 0.5 * np.einsum("mnef,efin->mi", Jmsa[o, o, v, v], ttau, optimize=True)
        Fme =                           Fms[o, v] +       np.einsum("mnef,fn->me",   Jmsa[o, o, v, v], t1,   optimize=True)

        # define some complementary variables used in the following expressions
        Fmea =            np.einsum("bm,me->be",   t1, Fme, optimize=True)
        Fmeb =            np.einsum("ej,me->mj",   t1, Fme, optimize=True)
        t12  = 0.5 * t2 + np.einsum("fj,bn->fbjn", t1, t1,  optimize=True)

        # define the permutation arguments for all terms the W intermediates
        P1 = np.einsum("ej,mnie->mnij", t1, Jmsa[o, o, o, v], optimize=True)
        P2 = np.einsum("bm,amef->abef", t1, Jmsa[v, o, v, v], optimize=True)

        # calculate the 4D two-particle intermediates
        Wmnij = Jmsa[o, o, o, o] + 0.25 * np.einsum("efij,mnef->mnij", tau, Jmsa[o, o, v, v], optimize=True) + P1 - P1.swapaxes(2, 3)
        Wabef = Jmsa[v, v, v, v] + 0.25 * np.einsum("abmn,mnef->abef", tau, Jmsa[o, o, v, v], optimize=True) - P2 + P2.swapaxes(0, 1)
        Wmbej = Jmsa[o, v, v, o] +        np.einsum("fj,mbef->mbej",   t1,  Jmsa[o, v, v, v], optimize=True)                                  -        np.einsum("bn,mnej->mbej",   t1,  Jmsa[o, o, v, o], optimize=True)                                  -        np.einsum("fbjn,mnef->mbej", t12, Jmsa[o, o, v, v], optimize=True)

        # define the right hand side of the T1 and T2 amplitude equations
        rhs_T1, rhs_T2 = Fms[v, o].copy(), Jmsa[v, v, o, o].copy()

        # calculate the right hand side of the CCSD equation for T1
        rhs_T1 +=       np.einsum("ei,ae->ai",     t1, Fae,              optimize=True)
        rhs_T1 -=       np.einsum("am,mi->ai",     t1, Fmi,              optimize=True)
        rhs_T1 +=       np.einsum("aeim,me->ai",   t2, Fme,              optimize=True)
        rhs_T1 -=       np.einsum("fn,naif->ai",   t1, Jmsa[o, v, o, v], optimize=True)
        rhs_T1 -= 0.5 * np.einsum("efim,maef->ai", t2, Jmsa[o, v, v, v], optimize=True)
        rhs_T1 -= 0.5 * np.einsum("aemn,nmei->ai", t2, Jmsa[o, o, v, o], optimize=True)

        # define the permutation arguments for all terms in the equation for T2
        P1  = np.einsum("aeij,be->abij",    t2,     Fae - 0.5 * Fmea, optimize=True)
        P2  = np.einsum("abim,mj->abij",    t2,     Fmi + 0.5 * Fmeb, optimize=True)
        P3  = np.einsum("aeim,mbej->abij",  t2,     Wmbej,            optimize=True)
        P3 -= np.einsum("ei,am,mbej->abij", t1, t1, Jmsa[o, v, v, o], optimize=True)
        P4  = np.einsum("ei,abej->abij",    t1,     Jmsa[v, v, v, o], optimize=True)
        P5  = np.einsum("am,mbij->abij",    t1,     Jmsa[o, v, o, o], optimize=True)

        # calculate the right hand side of the CCSD equation for T2
        rhs_T2 += 0.5 * np.einsum("abmn,mnij->abij", tau, Wmnij, optimize=True)
        rhs_T2 += 0.5 * np.einsum("efij,abef->abij", tau, Wabef, optimize=True)
        rhs_T2 += P1.transpose(0, 1, 2, 3) - P1.transpose(1, 0, 2, 3)
        rhs_T2 -= P2.transpose(0, 1, 2, 3) - P2.transpose(0, 1, 3, 2)
        rhs_T2 += P3.transpose(0, 1, 2, 3) - P3.transpose(0, 1, 3, 2)
        rhs_T2 -= P3.transpose(1, 0, 2, 3) - P3.transpose(1, 0, 3, 2)
        rhs_T2 += P4.transpose(0, 1, 2, 3) - P4.transpose(0, 1, 3, 2)
        rhs_T2 -= P5.transpose(0, 1, 2, 3) - P5.transpose(1, 0, 2, 3)

        # Update T1 and T2 amplitudes
        t1, t2 = rhs_T1 * Emss, rhs_T2 * Emsd

        # evaluate the energy
        E_CCSD_P, E_CCSD = E_CCSD, np.einsum("ia,ai", Fms[o, v], t1) + 0.25 * np.einsum("ijab,abij", Jmsa[o, o, v, v], t2) + 0.5 * np.einsum("ijab,ai,bj", Jmsa[o, o, v, v], t1, t1)

    # print the CCSD energy
    print("   CCSD ENERGY: {:.8f}".format(E_HF + E_CCSD + VNN))
```

{:.cite}
> Stanton, John F., Jürgen Gauss, John D. Watts, and Rodney J. Bartlett. 1991. “A Direct Product Decomposition Approach for Symmetry Exploitation in Many-Body Methods. I. Energy Calculations.” *The Journal of Chemical Physics* 94: 4334–45. <https://doi.org/10.1063/1.460620>.
>

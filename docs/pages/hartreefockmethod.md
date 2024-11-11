---
title: Hartree–Fock Method
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Hartree--Fock Method

The Hartree--Fock method is a foundational approach in quantum chemistry, aimed at solving the electronic structure problem in molecules by determining an optimal wavefunction. This method simplifies the complex interactions of electrons through a mean-field approximation, where each electron moves in an average field created by all others. This allows for the use of a single set of orbitals, leading to the construction of the Fock operator and iterative solutions to one-electron equations.

However, Hartree--Fock has notable limitations. Its reliance on a single-determinant wavefunction means it struggles to account for electron correlation. This leads to inaccuracies in energy predictions, particularly for systems with strong electron interactions, such as transition metal complexes or molecules with delocalized electrons.

## Theoretical Background

Our primary objective is to solve the Schrödinger equation in the form

\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi}
\end{equation}

where $\hat{\mathbf{H}}$ denotes the molecular Hamiltonian operator, $\ket{\Psi}$ represents the molecular wavefunction, and $E$ is the total energy of the system. The Hartree--Fock approximates the total wavefunction $\ket{\Psi}$ as a single Slater determinant, expressed as

\begin{equation}
\ket{\Psi}=\ket{\chi\_1\chi\_2\cdots\chi\_N}
\end{equation}

where $\chi\_i$ denotes a spin orbital, and $N$ is the total number of electrons. The goal of the Hartree--Fock method is to optimize these orbitals in order to minimize the system's total energy, thereby providing a reliable estimate of the electronic structure.

In the Restricted Hartree--Fock method, we impose a constraint on electron spin, which allows us to work with spatial orbitals instead of spin orbitals. This reformulation expresses the Slater determinant in terms of spatial orbitals as

\begin{equation}
\ket{\Psi}=\ket{\Phi\_1\Phi\_2\cdots\Phi\_{N/2}}
\end{equation}

where $\Phi\_i$ represents a spatial orbital. Notably, the Restricted Hartree--Fock method requires an even number of electrons to satisfy spin-pairing. In practice, atomic orbitals (whether spin or spatial) are typically expanded in a set of basis functions $\left\lbrace\phi\_i\right\rbrace$, which are often Gaussian functions, allowing for convenient computation with expansion coefficients. After performing some algebraic manipulations, the Hartree--Fock method leads to the Roothaan equations, which are expressed as

\begin{equation}\label{eq:roothaan}
\mathbf{FC}=\mathbf{SC}\mathbf{\varepsilon}
\end{equation}

where $\mathbf{F}$ is the Fock matrix, $\mathbf{C}$ is the matrix of orbital coefficients, $\mathbf{S}$ is the overlap matrix, and $\mathbf{\varepsilon}$ represents the orbital energies. These matrices will be defined in detail later.

## Implementation of the Restricted Hartree--Fock Method

To begin, we define the core Hamiltonian, also known as the one-electron Hamiltonian. This component of the full Hamiltonian excludes electron-electron repulsion and is expressed in index notation as

\begin{equation}\label{eq:hamiltonian}
H\_{\mu\nu}^{\mathrm{core}}=T\_{\mu\nu}+V\_{\mu\nu}
\end{equation}

where $\mu$ and $\nu$ are indices of the basis functions. Here, $T\_{\mu\nu}$ represents a kinetic energy matrix element, while $V\_{\mu\nu}$ denotes a potential energy matrix element. These matrix elements are defined as

\begin{align}
T\_{\mu\nu}&=\braket{\phi\_{\mu}|\hat{T}|\phi\_{\nu}} \\\\\
V\_{\mu\nu}&=\braket{\phi\_{\mu}|\hat{V}|\phi\_{\nu}}
\end{align}

Additionally, the overlap integrals, which describe the extent of overlap between basis functions, are given by

\begin{equation}\label{eq:overlap}
S\_{\mu\nu}=\braket{\phi\_{\mu}|\phi\_{\nu}}
\end{equation}

and another essential component are the two-electron repulsion integrals, defined as

\begin{equation}\label{eq:coulomb}
J\_{\mu\nu\kappa\lambda}=\braket{\phi\_{\mu}\phi\_{\mu}|\hat{J}|\phi\_{\kappa}\phi\_{\lambda}}
\end{equation}

which play crucial roles in the Hartree--Fock calculation. All of these integrals over (Gausssian) basis functions are usually calculated using analytical expressions.<!--\cite{10.1016/S0065-3276!08!60019-2}--> The solution to the Roothaan equations \eqref{eq:roothaan} requires an iterative procedure, since the Fock matrix defined as

\begin{equation}\label{eq:fock}
F\_{\mu\nu}=H\_{\mu\nu}^{\mathrm{core}}+D\_{\kappa\lambda}\left(J\_{\mu\nu\kappa\lambda}-\frac{1}{2}J\_{\mu\lambda\kappa\nu}\right)
\end{equation}

depends on the unknown density matrix $\mathbf{D}$, defined as

\begin{equation}\label{eq:hf_density}
D\_{\mu\nu}=2C\_{\mu i}C\_{\nu i}
\end{equation}

This iterative procedure is executed through the Self-Consistent Field method. An initial guess is made for the density matrix $\mathbf{D}$ (usually zero), the Roothaan equations \eqref{eq:roothaan} are solved and the density matix is updated using the equation \eqref{eq:hf_density}. The total energy of the system is then calculated using the core Hamiltonian and the Fock matrix as

\begin{equation}
E=\frac{1}{2}D\_{\mu\nu}(H\_{\mu\nu}^{\mathrm{core}}+F\_{\mu\nu})+E\_{\mathrm{nuc}}
\end{equation}

After convergence of both the density matrix and total energy, the process concludes, yielding the optimized molecular orbitals. The total energy of the system also includes the nuclear repulsion energy, which is given by

\begin{equation}
E\_{\mathrm{nuc}}=\sum\_{A}\sum\_{B<A}\frac{Z\_{A}Z\_{B}}{R\_{AB}}
\end{equation}

where $Z\_A$ is the nuclear charge of atom $A$, and $R\_{AB}$ is the distance between atoms $A$ and $B$.

### Direct Inversion in the Iterative Subspace

In the Hartree--Fock method, convergence of the density matrix and energy can be significantly accelerated by employing the Direct Inversion in the Iterative Subspace technique. Direct Inversion in the Iterative Subspace achieves this by storing Fock matrices from previous iterations and constructing an optimized linear combination that minimizes the current iteration's error. This approach is especially valuable in Hartree--Fock calculations for large systems, where convergence issues are more common and challenging to resolve.

We start by defining the error vector of $i$-th iteration $\mathbf{e}\_i$ as

\begin{equation}
\mathbf{e}\_i=\mathbf{S}\_i\mathbf{D}\_i\mathbf{F}\_i-\mathbf{F}\_i\mathbf{D}\_i\mathbf{S}\_i
\end{equation}

Our goal is to transform the Fock matrix $\mathbf{F}\_i$ as

\begin{equation}\label{eq:fock_extrapolate}
\mathbf{F}\_i=\sum\_{j=i-(L+1)}^{i-1}c\_j\mathbf{F}\_j
\end{equation}

where $c\_j$ are the coefficients that minimize the error matrix $\mathbf{e}\_i$ and $L$ is the number of Fock matrices we store called the subspace size. To calculate the coefficients $c\_j$, we solve the set linear equations

\begin{equation}
\begin{bmatrix}
\mathbf{e}\_1\cdot\mathbf{e}\_1 & \dots & \mathbf{e}\_1\cdot\mathbf{e}\_{L+1} & 1 \\\\\
\vdots & \ddots & \vdots & \vdots \\\\\
\mathbf{e}\_{L+1}\cdot\mathbf{e}\_1 & \dots & \mathbf{e}\_{L+1}\cdot\mathbf{e}\_{L+1} & 1 \\\\\
1 & \dots & 1 & 0 \\\\\
\end{bmatrix}
\begin{bmatrix}
c\_1 \\\\\
\vdots \\\\\
c\_{L+1} \\\\\
\lambda \\\\\
\end{bmatrix}
=
\begin{bmatrix}
0 \\\\\
\vdots \\\\\
0 \\\\\
1 \\\\\
\end{bmatrix}
\end{equation}

After solving the linear equations, we use the coefficients $c\_j$ to construct the new Fock matrix $\mathbf{F}\_i$ according to the equation \eqref{eq:fock_extrapolate} and proceed as usual with the Hartree--Fock calculation.

### Gradient of the Restricted Hartree--Fock Method

If we perform the calculation as described above and get the density matrix $\mathbf{D}$ we can evaluate the nuclear energy gradient as<!--\cite{10.1002/9780470749593.hrs006}-->

\begin{equation}\label{eq:hf_gradient}
\frac{\partial E}{\partial X\_{A,i}}=\sum\_{\mu\nu\in\left\lbrace\phi\_A\right\rbrace}D\_{\mu\nu}\frac{\partial H\_{\mu\nu}^{\mathrm{core}}}{\partial X\_{A,i}}+2\sum\_{\mu\nu\in\left\lbrace\phi\_A\right\rbrace}D\_{\mu\nu}D\_{\kappa\lambda}\frac{\partial\left(J\_{\mu\nu\kappa\lambda}-\frac{1}{2}J\_{\mu\lambda\kappa\nu}\right)}{\partial X\_{A,i}}-2W\_{\mu\nu}\frac{\partial S\_{\mu\nu}}{\partial X\_{A,i}}
\end{equation}

where $A$ is the index of an atom, $i$ is the index of the coordinate, $\left\lbrace\phi\_A\right\rbrace$ is the set of all basis functions located at atom $A$ and $\mathbf{W}$ is energy weighed density matrix defined as

\begin{equation}
W_{\mu\nu}=2C_{\mu i}C_{\nu i}\varepsilon_i
\end{equation}

Keep in mind that the indices $\kappa$ and $\lambda$ in the gradient equation \eqref{eq:hf_gradient} are summed over all basis functions.

## Integral Transforms to the Basis of Molecular Spinorbitals<!--\label{sec:integral_transform}-->

To carry out most post-Hartree--Fock calculations, it is essential to transform the integrals into the Molecular Spinorbital basis. We will outline this transformation process here and refer to it in subsequent sections on post-Hartree--Fock methods. The post-Hartree--Fock methods in this document will be presented using the integrals in the Molecular Spinorbital basis (and in its antisymmetrized form for the two-electron integrals) as this approach is more general.

All the integrals defined in the equations \eqref{eq:hamiltonian}, \eqref{eq:overlap}, and \eqref{eq:coulomb}, as well as Fock matrix in the equation \eqref{eq:fock} are defined in the basis of atomic orbitals. To transform these integrals to the Molecular Spinorbital basis, we utilize the coefficient matrix $\mathbf{C}$ obtained from the solution of the Roothaan equations \eqref{eq:roothaan}. This coefficient matrix $\mathbf{C}$ is initially calculated in the spatial molecular orbital basis (in the Restricted Hartree--Fock calculation).

The first step involves expanding the coefficient matrix $\mathbf{C}$ to the Molecular Spinorbital basis. This transformation can be mathematically expressed using the tiling matrix $\mathbf{P}\_{n\times 2n}$, defined as

\begin{equation}
\mathbf{P}=
\begin{pmatrix}
e\_1&e\_1&e\_2&e\_2&\dots&e\_n&e\_n
\end{pmatrix}
\end{equation}

where $e\_i$ represents the $i$-th column of the identity matrix $\mathbf{I}\_n$. Additionally, we define the matrices $\mathbf{M}\_{n\times 2n}$ and $\mathbf{N}\_{n\times 2n}$ with elements given by

\begin{equation}
M\_{ij}=1-j\bmod 2,N\_{ij}=j \bmod 2
\end{equation}

The coefficient matrix $\mathbf{C}$ in the Molecular Spinorbital basis can be then expressed as

\begin{equation}
\mathbf{C}^{\mathrm{MS}}=
\begin{pmatrix}
\mathbf{CP} \\\\\
\mathbf{CP}
\end{pmatrix}
\odot
\begin{pmatrix}
\mathbf{M} \\\\\
\mathbf{N}
\end{pmatrix}
\end{equation}

where $\odot$ denotes the Hadamard product. This transformed matrix $\mathbf{C}^{\mathrm{MS}}$ is subsequently used to transform the two-electron integrals $\mathbf{J}$ to the Molecular Spinorbital basis as

\begin{equation}
J\_{pqrs}^{\mathrm{MS}}=C\_{\mu p}^{\mathrm{MS}}C\_{\nu q}^{\mathrm{MS}}(\mathbf{I}\_{2}\otimes\_K(\mathbf{I}\_{2}\otimes\_K\mathbf{J})^{(4,3,2,1)})\_{\mu\nu\kappa\lambda}C\_{\kappa r}^{\mathrm{MS}}C\_{\lambda s}^{\mathrm{MS}}
\end{equation}

where the superscript $(4,3,2,1)$ denotes the axes transposition and $\otimes\_K$ is the Kronecker product. This notation accommodates the spin modifications and ensures adherence to quantum mechanical principles. We also define the antisymmetrized two-electron integrals in physicists' notation as

\begin{equation}
\braket{pq||rs}=(J\_{pqrs}^{\mathrm{MS}}-J\_{psrq}^{\mathrm{MS}})^{(1,3,2,4)}
\end{equation}

For the transformation of the one-electron integrals such as the core Hamiltonian, the overlap matrix and also the Fock matrix, we use the formula

\begin{equation}
A\_{pq}^{\mathrm{MS}}=C\_{\mu p}^{\mathrm{MS}}(\mathbf{I}\_{2}\otimes\_K\mathbf{A})\_{\mu\nu}C\_{\nu q}^{\mathrm{MS}}
\end{equation}

where $\mathbf{A}$ is an arbitrary matrix of one-electron integrals. Since many post-Hartree--Fock methods rely on differences of orbital energies in the denominator, we define the tensors

\begin{align}
\varepsilon^{a}\_{i}&=\varepsilon\_i-\varepsilon\_a \\\\\
\varepsilon^{ab}\_{ij}&=\varepsilon\_i+\varepsilon\_j-\varepsilon\_a-\varepsilon\_b \\\\\
\varepsilon^{abc}\_{ijk}&=\varepsilon\_i+\varepsilon\_j+\varepsilon\_k-\varepsilon\_a-\varepsilon\_b-\varepsilon\_c
\end{align}

for convenience. These tensors enhance code readability and efficiency, making it easier to understand and work with the underlying mathematical framework. Here and also throughout the rest of the document, the indices $i$, $j$ and $k$ run over occupied orbitals, whereas the indices $a$, $b$ and $c$ run over virtual orbitals.

## Hartree--Fock Method and Integral Transform Coding Exercise<!--\label{sec:hf_int_code_exercise}-->

This section provides Python code snippets for implementing the Hartree--Fock method and transforming integrals to the Molecular Spinorbital basis. The code uses the NumPy library for efficient numerical computations, and the exercises are designed to build familiarity with the Hartree--Fock method and the Molecular Spinorbital integral transformations.

Each exercise includes placeholders for you to fill. It is assumed that the variables `atoms`, `coords`, `S`, `H`, and `J` are defined before the exercises. These variables represent the atomic numbers, atomic coordinates, overlap matrix, core Hamiltonian, and two-electron integrals, respectively. These variables can typically be obtained from the output of a quantum chemistry software package. If you would like to focus solely on coding you can save the [molecule file](/acorn/python/molecule.xyz), [overlap integrals](/acorn/python/S_AO.mat), [core Hamiltonian](/acorn/python/H_AO.mat), and [two-electron integrals](/acorn/python/J_AO.mat) in the STO-3G basis to the same directory as the exercise code and load the variables using the Listing <!--\ref{code:load_exercise}--> below. The `ATOM` variable is a dictionary that maps the atomic symbols to atomic numbers.

<!--{id=code:load_exercise caption="Example loading of molecule and integrals over atomic basis functions into variables used throughout exercises. The snippet expects the molecule and integral files to b present in the same directory as the script."}-->
```py
# get the atomic numbers and coordinates of all atoms
atoms = np.array([ATOM[line.split()[0]] for line in open("molecule.xyz").readlines()[2:]], dtype=int)
coords = np.array([line.split()[1:] for line in open("molecule.xyz").readlines()[2:]], dtype=float)

# convert to bohrs
coords *= 1.8897261254578281

# load the integrals from the files
H, S = np.loadtxt("H_AO.mat", skiprows=1), np.loadtxt("S_AO.mat", skiprows=1); J = np.loadtxt("J_AO.mat", skiprows=1).reshape(4 * [S.shape[1]])
```

With all the variables defined, you can proceed to the Hartree--Fock and integral transform exercises in the Listings <!--\ref{code:hf_exercise} and \ref{code:int_exercise}--> below.

<!--{id=code:hf_exercise caption="Hartree--Fock method exercise code. The important variables like number of occupied orbitals, convergence threshold and matrix containers are already defined. The student is expected to fill the Self-Consistent Field loop and calculate nuclear repulsion energy from the atomic numbers and coordinates. The total energy is then automatically printed."}-->
```python
"""
Here are defined some of the necessary variables. The variable "E_HF" stores the Hartree-Fock energy, while "E_HF_P" keeps track of the previous iteration's energy to monitor convergence. The "thresh" defines the convergence criteria for the calculation. The variables "nocc" and "nbf" represent the number of occupied orbitals and the number of basis functions, respectively. Initially, "E_HF" is set to zero and "E_HF_P" to one to trigger the start of the Self-Consistent Field (SCF) loop. Although you can rename these variables, it is important to note that certain sections of the code are tailored to these specific names.
"""
E_HF, E_HF_P, nocc, nbf, thresh = 0, 1, sum(atoms) // 2, S.shape[0], 1e-8

"""
These lines set up key components for our HF calculations. We initialize the density matrix as a zero matrix, and the coefficients start as an empty array. Although the coefficient matrix is computed within the while loop, it's defined outside to allow for its use in subsequent calculations, such as the MP energy computation. Similarly, the exchange tensor is accurately calculated here by transposing the Coulomb tensor. The "eps" vector, which contains the orbital energies, is also defined at this stage to facilitate access throughout the script. This setup ensures that all necessary variables are ready for iterative processing and further calculations beyond the SCF loop.
"""
K, F, D, C, eps = J.transpose(0, 3, 2, 1), np.zeros_like(S), np.zeros_like(S), np.zeros_like(S), np.array(nbf * [0])

"""
This while loop is the SCF loop. Please fill it so it calculates the Fock matrix, solves the Fock equations, builds the density matrix from the coefficients and calculates the energy. You can use all the variables defined above and all the functions in numpy package. The recommended functions are np.einsum and np.linalg.eigh. Part of the calculation will probably be calculation of the inverse square root of a matrix. The numpy package does not conatin a function for this. You can find a library that can do that or you can do it manually. The manual calculation is, of course, preferred.
"""
while abs(E_HF - E_HF_P) > thresh:
    break

"""
In the followng block of code, please calculate the nuclear-nuclear repulsion energy. You should use only the atoms and coords variables. The code can be as short as two lines. The result should be stored in the "VNN" variable.
"""
VNN = 0

# print the results
print("RHF ENERGY: {:.8f}".format(E_HF + VNN))
```

<!--{id=code:int_exercise caption="Integral transform exercise code. The tiling matrices are predefined here with a correct shape, but the student is expected to fill the with the correct expressions. After that, the student should transform the coefficient matrix, core Hamoltonian, Fock matrix and two-electron integrals to the Molecular Spinorbital basis. Additionally, orbital energy tensor are also expected to be calculated."}-->
```python
"""
To perform most of the post-HF calculations, we need to transform the Coulomb integrals to the molecular spinorbital basis, so if you don't plan to calculate any post-HF methods, you can end the eercise here. The restricted MP2 calculation could be done using the Coulomb integral in MO basis, but for the sake of subsequent calculations, we enforce here the integrals in the MS basis. The first thing you will need for the transform is the coefficient matrix in the molecular spinorbital basis. To perform this transform using the mathematical formulation presented in the materials, the first step is to form the tiling matrix "P" which will be used to duplicate columns of a general matrix. Please define it here.
"""
P = np.zeros((nbf, 2 * nbf))

"""
Now, please define the spin masks "M" and "N". These masks will be used to zero out spinorbitals, that should be empty.
"""
M, N = np.zeros((nbf, 2 * nbf)), np.zeros((nbf, 2 * nbf))

"""
With the tiling matrix and spin masks defined, please transform the coefficient matrix into the molecular spinorbital basis. The resulting matrix should be stored in the "Cms" variable.
"""
Cms = np.zeros(2 * np.array(C.shape))

"""
For some of the post-HF calculations, we will also need the Hamiltonian and Fock matrix in the molecular spinorbital basis. Please transform it and store it in the "Hms" and "Fms" variable. If you don't plan to calculate the CCSD method, you can skip the transformation of the Fock matrix, as it is not needed for the MP2 and CI calculations.
"""
Hms, Fms = np.zeros(2 * np.array(H.shape)), np.zeros(2 * np.array(H.shape))

"""
With the coefficient matrix in the molecular spinorbital basis available, we can proceed to transform the Coulomb integrals. It is important to note that the transformed integrals will contain twice as many elements along each axis compared to their counterparts in the atomic orbital (AO) basis. This increase is due to the representation of both spin states in the molecular spinorbital basis.
"""
Jms = np.zeros(2 * np.array(J.shape))

"""
The post-HF calculations also require the antisymmetrized two-electron integrals in the molecular spinorbital basis. These integrals are essential for the MP2 and CC calculations. Please define the "Jmsa" tensor as the antisymmetrized two-electron integrals in the molecular spinorbital basis.
"""
Jmsa = np.zeros(2 * np.array(J.shape))

"""
As mentioned in the materials, it is also practical to define the tensors of reciprocal orbital energy differences in the molecular spinorbital basis. These tensors are essential for the MP2 and CC calculations. Please define the "Emss", "Emsd" and "Emst" tensors as tensors of single, double and triple excitation energies, respectively. The configuration interaction will not need these tensors, so you can skip this step if you don't plan to program the CI method. The MP methods will require only the "Emsd" tensor, while the CC method will need both tensors.
"""
Emss, Emsd = np.array([]), np.array([])
```

Solution to this exercises can be found [here](codesolutions.html#hartreefock-method) and [here](codesolutions.html#integral-transform).

{:.cite}
> Gill, Peter M. W. 1994. *Molecular Integrals over Gaussian Basis Functions*. Academic Press. <https://doi.org/10.1016/S0065-3276(08)60019-2>.
>

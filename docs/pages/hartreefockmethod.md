---
title: Hartree–Fock Method
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Hartree--Fock Method

The Hartree--Fock method stands as a cornerstone in quantum chemistry, offering a systematic approach to solve the electronic structure problem in molecules. This computational technique strives to determine the optimal wave function for a given molecular system, providing insights into the distribution of electrons and their energies.

## Theoretical Background

Ultimately, we are interested in solving the Schrödinger equation in the form

\begin{equation}
\hat{\mathbf{H}}\ket{\Psi}=E\ket{\Psi},
\end{equation}

where $\hat{H}$ is the molecular Hamiltonian operator, $\ket{\Psi}$ is the molecular wave function, and $E$ is the total energy of the system. The Hartree--Fock method aims to approximate the wave function $\ket{\Psi}$ by a single Slater determinant, which we will write in the form

\begin{equation}
\ket{\Psi}=\ket{\chi\_1\chi\_2\cdots\chi\_N},
\end{equation}

where $\chi\_i$ represents a Molecular Spinorbital and $N$ is the total number of electrons. The Hartree--Fock method seeks to optimize these molecular orbitals to minimize the total energy of the system, providing a reliable estimate of the electronic structure. However, in the Restricted Hartree--Fock method, we impose a constraint on the electron spin, and instead of spin orbitals, we work with spatial orbitals, which allows us to write the slater determinant in terms of spatial orbitals in the form

\begin{equation}
\ket{\Psi}=\ket{\Phi\_1\Phi\_2\cdots\Phi\_{N/2}},
\end{equation}

where $\Phi\_i$ represents a molecular spatial orbital. We can see, that for the Restricted Hartree--Fock method, we need to have an even number of electrons in the system. In practical calculations, it is convenient to expand the molecular orbitals (spin or spatial) in terms of basis functions $\phi$ (usually sums of Gaussian functions) and work with the expansion coefficients. If we assume that the wavefunction is a single Slater determinant, our molecular orbitals are expanded in terms of basis functions and optimize the energy of such determinant, we arrive at the Roothaan equations in the form

\begin{equation}\label{eq:roothaan}
\mathbf{FC}=\mathbf{SC\varepsilon},
\end{equation}

where $\mathbf{F}$ is the Fock matrix (defined later), $\mathbf{C}$ is a matrix of orbital coefficients, $\mathbf{S}$ is the overlap matrix (also defined later), and $\mathbf{\varepsilon}$ represents orbital energies.

## Implementation of the Restricted Hartree–Fock Method

Let's begin by defining the core Hamiltonian, also known as the one-electron Hamiltonian. The core Hamiltonian represents a part of the full Hamiltonian that excludes electron-electron repulsion. In index notation, it is expressed as 

\begin{equation}\label{eq:hamiltonian}
H\_{\mu\nu}^{core}=T\_{\mu\nu}+V\_{\mu\nu}
\end{equation}

where $\mu$ and $\nu$ are indices of basis functions, $T\_{\mu\nu}$ is a kinetic energy matrix element and $V\_{\mu\nu}$ is a potential energy matrix element. These matrix elements are given as

\begin{align}
T\_{\mu\nu}&=\braket{\phi\_{\mu}|\hat{T}|\phi\_{\nu}} \\\\\
V\_{\mu\nu}&=\braket{\phi\_{\mu}|\hat{V}|\phi\_{\nu}}
\end{align}

are usually calculated using analytical expressions. Additionally, using analytical expressions, we can calculate the overlap integrals

\begin{equation}\label{eq:overlap}
S\_{\mu\nu}=\braket{\phi\_{\mu}|\phi\_{\nu}}
\end{equation}

and the two-electron Coulomb repulsion integrals

\begin{equation}\label{eq:coulomb}
J\_{\mu\nu\kappa\lambda}=\braket{\phi\_{\mu}\phi\_{\mu}|\hat{J}|\phi\_{\kappa}\phi\_{\lambda}},
\end{equation}

which play crucial roles in the Hartree--Fock calculation.<!--\cite{10.1016/S0065-3276!08!60019-2}--> The Hartree--Fock method revolves around solving the Roothaan equations \eqref{eq:roothaan} iteratively. The Fock matrix is defined as

\begin{equation}\label{eq:fock}
F\_{\mu\nu}=H\_{\mu\nu}^{core}+D\_{\kappa\lambda}(J\_{\mu\nu\kappa\lambda}-\frac{1}{2}J\_{\mu\lambda\kappa\nu})
\end{equation}

depends on the unknown density matrix $\mathbf{D}$. This iterative process is carried out through a Self-Consistent Field method. In each iteration, a guess for the density matrix is made, and the Roothaan equations are solved. The density matrix is then updated as

\begin{equation}
D\_{\mu\nu}=2C\_{\mu i}C\_{\nu i}
\end{equation}

and the total energy of the system

\begin{equation}
E=\frac{1}{2}D\_{\mu\nu}(H\_{\mu\nu}^{core}+F\_{\mu\nu})+E\_{nuc}
\end{equation}

is then calculated using the core Hamiltonian and the Fock matrix. The important thing to note is that the Fock matrix depends on the density matrix, which is updated in each iteration. The initial guess for the density matrix is often set to zero. After the density matrix and the total energy are converged, the Self-Consistent Field procedure is terminated, and the optimized molecular orbitals are obtained. To get the final energy of the system, we need to add the nuclear repulsion energy of the form

\begin{equation}
E\_{nuc}=\sum\_{A}\sum\_{B<A}\frac{Z\_{A}Z\_{B}}{R\_{AB}},
\end{equation}

where $Z\_A$ is the nuclear charge of atom A, and $R\_{AB}$ is the distance between atoms A and B.

### Gradient of the Restricted Hartree--Fock Method

If we perform the calculation as described above and get the density matrix $\mathbf{D}$ we can evaluate the nuclear energy gradient as

\begin{equation}
\frac{\partial E}{\partial X\_{A,i}}=D\_{\mu\nu}\frac{\partial H\_{\mu\nu}^{core}}{\partial X\_{A,i}}+2D\_{\mu\nu}D\_{\kappa\lambda}\frac{\partial J\_{\mu\nu\kappa\lambda}}{\partial X\_{A,i}}-2W\_{\mu\nu}\frac{\partial S\_{\mu\nu}}{\partial X\_{A,i}}
\end{equation}

where $i$ is the index of the coordinate and where $\mathbf{W}$ is energy weighed density matrix defined as

\begin{equation}
W_{\mu\nu}=2C_{\mu i}C_{\nu i}\varepsilon_i
\end{equation}

## Integral Transforms to the Basis of Molecular Spinorbitals<!--\label{sec:integral_transform}-->

To perform most of the post-Hartree--Fock calculations, we usually need to transform the integrals to the Molecular Spinorbital basis. We will describe it here and refer to it in the post-Hartree--Fock methods sections. We will also present the post-Hartree--Fock methods using the integrals in the Molecular Spinorbital basis (and its antisymmetrized form in case of the Coulomb integrals), since it is more general.

All the integrals defined in the equations \eqref{eq:hamiltonian}, \eqref{eq:overlap}, and \eqref{eq:coulomb} and even the Fock matrix in the equation \eqref{eq:fock} are defined in the basis of atomic orbitals. To transform these integrals to the Molecular Spinorbital basis, we need to use the coefficient matrix $\mathbf{C}$ obtained from the solution of the Roothaan equations \eqref{eq:roothaan}. The coefficient matrix $\mathbf{C}$, which is obtained from the Restricted Hartree--Fock calculation, is calculated in the spatial molecular orbital basis. The first step is to expand the coefficient matrix $\mathbf{C}$ to the Molecular Spinorbital basis. This can be done mathematically using the tiling matrix $\mathbf{P}\_{n\times 2n}$, defined as

\begin{equation}
\mathbf{P}=
\begin{pmatrix}
e\_1&e\_1&e\_2&e\_2&\dots&e\_n&e\_n
\end{pmatrix}
,
\end{equation}

where $e\_i$ represents the $i$-th column of the identity matrix $\mathbf{I}\_n$ and the matrices $\mathbf{M}\_{n\times 2n}$ and $\mathbf{N}\_{n\times 2n}$ with elements given by

\begin{equation}
M\_{ij}=1-j\bmod 2,N\_{ij}=j \bmod 2.
\end{equation}

The coefficient matrix $\mathbf{C}$ in the MS basis can be then expressed as

\begin{equation}
\mathbf{C}^{MS}=
\begin{pmatrix}
\mathbf{CP} \\\\\
\mathbf{CP}
\end{pmatrix}
\odot
\begin{pmatrix}
\mathbf{M} \\\\\
\mathbf{N}
\end{pmatrix}
,
\end{equation}

where $\odot$ denotes the Hadamard product. This transformed matrix $\mathbf{C}^{MS}$ is then used to transform the Coulomb integrals $\mathbf{J}$ to the MS basis as

\begin{equation}
J\_{pqrs}^{MS}=C\_{\mu p}^{MS}C\_{\nu q}^{MS}(\mathbf{I}\_{2}\otimes\_K(\mathbf{I}\_{2}\otimes\_K\mathbf{J})^{(4,3,2,1)})\_{\mu\nu\kappa\lambda}C\_{\kappa r}^{MS}C\_{\lambda s}^{MS},
\end{equation}

where the superscript $(4,3,2,1)$ denotes the axes transposition and $\otimes\_K$ is the Kronecker product. This notation accounts for the spin modifications and ensures that the transformations adhere to quantum mechanical principles. We also define the antisymmetrized Coulomb integrals in physicists' notation as

\begin{equation}
\braket{pq||rs}=(J\_{pqrs}^{MS}-J\_{psrq}^{MS})^{(1,3,2,4)}.
\end{equation}

For the transformation of the one-electron integrals such as the core Hamiltonian, the overlap matrix and also the Fock matrix, we use the formula

\begin{equation}
A\_{pq}^{MS}=C\_{\mu p}^{MS}(\mathbf{I}\_{2}\otimes\_K\mathbf{A})\_{\mu\nu}C\_{\nu q}^{MS},
\end{equation}

where $\mathbf{A}$ is an arbitrary one-electron integral. Since a lot of the post-Hartree--Fock methods also use differences of orbital energies in the denominator, it is practical to define the tensors

\begin{align}
\Delta\varepsilon^{a}\_{i}&=\frac{1}{\varepsilon\_i-\varepsilon\_a} \\\\\
\Delta\varepsilon^{ab}\_{ij}&=\frac{1}{\varepsilon\_i+\varepsilon\_j-\varepsilon\_a-\varepsilon\_b} \\\\\
\Delta\varepsilon^{abc}\_{ijk}&=\frac{1}{\varepsilon\_i+\varepsilon\_j+\varepsilon\_k-\varepsilon\_a-\varepsilon\_b-\varepsilon\_c},
\end{align}

where $a$, $b$, $c$ are virtual orbitals and $i$, $j$, $k$ are occupied orbitals. These tensors make the code more readable, easier to understand and also more efficient.

## Hartree--Fock Method and Integral Transform Code Example<!--\label{sec:hf_code_examples}-->

This section provides code snippets for the Hartree--Fock method and the integral transforms to the Molecular Spinorbital basis. The code snippets are written in Python and use the NumPy package for numerical calculations. The exercises are designed to help you understand the implementation of the Hartree--Fock method and the transformation of integrals to the Molecular Spinorbital basis. The solutions are provided to guide you through the implementation process and ensure that you can verify your results.

### Exercise

The exercise codes are designed to be self-contained and can be run in any Python environment. They contain placeholders that you need to fill in to complete the implementation. The exercises assume that you have defined the `atoms`, `coords`, `S`, `H`, and `J` variables, which represent the list of atomic numbers, atomic coordinates, overlap matrix, core Hamiltonian, and Coulomb integral tensor, respectively. These variables can be generaly obtained from a .xyz molecule file and the output of a quantum chemistry software. If you want to just get to the coding part, you can save the [molecule.xyz](/acorn/python/molecule.xyz), [S_AO.mat](/acorn/python/S_AO.mat), [H_AO.mat](/acorn/python/H_AO.mat), and [J_AO.mat](/acorn/python/J_AO.mat) files to the same directory as the exercise codes and and load the variables using the Listing below. The `ATOM` variable is a dictionary that maps the atomic symbols to atomic numbers.

```py
# get the atomic numbers and coordinates of all atoms
atoms = np.array([ATOM[line.split()[0]] for line in open("molecule.xyz").readlines()[2:]], dtype=int)
coords = np.array([line.split()[1:] for line in open("molecule.xyz").readlines()[2:]], dtype=float)

# convert to bohrs
coords *= 1.8897261254578281

# load the integrals from the files
H, S = np.loadtxt("H_AO.mat", skiprows=1), np.loadtxt("S_AO.mat", skiprows=1); J = np.loadtxt("J_AO.mat", skiprows=1).reshape(4 * [S.shape[1]])
```

With all the variables defined, you can proceed to the Hartree--Fock exercise in the Listing <!--\ref{code:hf_exercise}--> below.

<!--{id=code:hf_exercise caption="Hartree--Fock method exercise code."}-->
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

Now is a good time to check your results with the provided solution. If you are satisfied with the results, you can proceed to the next exercise in the Listing <!--\ref{code:int_exercise}--> below, which involves transforming the integrals to the molecular spinorbital basis.

<!--{id=code:int_exercise caption="Integral transform exercise code."}-->
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

If you successfully completed the exercise, you can compare your results with the provided solution. If you are satisfied with the results, you are now set for the post-HF methods exercises.

### Solution

The solutions provided below are complete implementations of the Hartree--Fock method <!--(Listing \ref{code:hf_solution})--> and the integral transforms <!--(Listing \ref{code:int_solution})--> to the Molecular Spinorbital basis. The solutions are written in Python and use the NumPy package for numerical calculations. The solutions are designed to guide you through the implementation process and ensure that you can verify your results.

<!--{id=code:hf_solution caption="Hartree--Fock method exercise code solution."}-->
```python
# define energies, number of occupied and virtual orbitals and the number of basis functions
E_HF, E_HF_P, VNN, nocc, nvirt, nbf = 0, 1, 0, sum(atoms) // 2, S.shape[0] - sum(atoms) // 2, S.shape[0]

# define some matrices and tensors
K, F, D, C, eps = J.transpose(0, 3, 2, 1), np.zeros_like(S), np.zeros_like(S), np.zeros_like(S), np.array(nbf * [0])

# calculate the X matrix which is the inverse of the square root of the overlap matrix
SEP = np.linalg.eigh(S); X = SEP[1] @ np.diag(1 / np.sqrt(SEP[0])) @ SEP[1].T

# the scf loop
while abs(E_HF - E_HF_P) > args.threshold:

    # build the Fock matrix
    F = H + np.einsum("ijkl,ij->kl", J - 0.5 * K, D, optimize=True)

    # solve the Fock equations
    eps, C = np.linalg.eigh(X @ F @ X); C = X @ C

    # build the density from coefficients
    D = 2 * np.einsum("ij,kj->ik", C[:, :nocc], C[:, :nocc])

    # save the previous energy and calculate the current electron energy
    E_HF_P, E_HF = E_HF, 0.5 * np.einsum("ij,ij", D, H + F, optimize=True)

# calculate nuclear-nuclear repulsion
for i, j in it.product(range(len(atoms)), range(len(atoms))):
    VNN += 0.5 * atoms[i] * atoms[j] / np.linalg.norm(coords[i, :] - coords[j, :]) if i != j else 0

# print the results
print("    RHF ENERGY: {:.8f}".format(E_HF + VNN))
```

<!--{id=code:int_solution caption="Integral transform exercise code solution."}-->
```python
# define the occ and virt spinorbital slices shorthand
o, v = slice(0, 2 * nocc), slice(2 * nocc, 2 * nbf)

# define the tiling matrix for the Jmsa coefficients and energy placeholders
P = np.array([np.eye(nbf)[:, i // 2] for i in range(2 * nbf)]).T

# define the spin masks
M = np.repeat([1 - np.arange(2 * nbf, dtype=int) % 2], nbf, axis=0)
N = np.repeat([    np.arange(2 * nbf, dtype=int) % 2], nbf, axis=0)

# tile the coefficient matrix, apply the spin mask and tile the orbital energies
Cms, epsms = np.block([[C @ P], [C @ P]]) * np.block([[M], [N]]), np.repeat(eps, 2)

# transform the core Hamiltonian and Fock matrix to the molecular spinorbital basis
Hms = np.einsum("ip,ij,jq->pq", Cms, np.kron(np.eye(2), H), Cms, optimize=True)
Fms = np.einsum("ip,ij,jq->pq", Cms, np.kron(np.eye(2), F), Cms, optimize=True)

# transform the coulomb integrals to the MS basis in chemist's notation
Jms = np.einsum("ip,jq,ijkl,kr,ls->pqrs", Cms, Cms, np.kron(np.eye(2), np.kron(np.eye(2), J).T), Cms, Cms, optimize=True);

# antisymmetrized two-electron integrals in physicist's notation
Jmsa = (Jms - Jms.swapaxes(1, 3)).transpose(0, 2, 1, 3)

# create the tensors of reciprocal differences of orbital energies in MS basis used in post-HF methods
Emss, Emsd = 1 / (epsms[o].reshape(-1) - epsms[v].reshape(-1, 1)), 1 / (epsms[o].reshape(-1) + epsms[o].reshape(-1, 1) - epsms[v].reshape(-1, 1, 1) - epsms[v].reshape(-1, 1, 1, 1))
```

{:.cite}
> Gill, Peter M. W. 1994. *Molecular Integrals over Gaussian Basis Functions*. Academic Press. <https://doi.org/10.1016/S0065-3276(08)60019-2>.
>

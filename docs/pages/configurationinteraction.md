---
title: Configuration Interaction
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Configuration Interaction

Configuration interaction (CI) is a post-Hartree-Fock, utilizing a linear variational approach to address the nonrelativistic Schrödinger equation under the Born–Oppenheimer approximation for multi-electron quantum systems. CI mathematically represents the wave function as a linear combination of Slater determinants. The term "configuration" refers to different ways electrons can occupy orbitals, while "interaction" denotes the mixing of these electronic configurations or states. CI computations, however, are resource-intensive, requiring significant CPU time and memory, limiting their application to smaller molecular systems. While Full CI (FCI) considers all possible electronic configurations, making it computationally prohibitive for larger systems, truncated versions like CISD, CISDT, and CISDTQ are more feasible and commonly employed in quantum chemistry studies.


## Restricted Full Configuration Interaction

Let's consider the restricted version of the Full Configuration Interaction (FCI) method, which considers all possible electronic configurations within a given basis set. The FCI wave function is expressed as a linear combination of Slater determinants, where each determinant represents a unique electron configuration. The FCI method provides the most accurate description of the electronic structure, but its computational cost grows exponentially with the number of electrons and basis functions, making it infeasible for large systems.

### The Transformation of Integrals to Molecular Spinorbital Basis
To begin, we need to convert the coefficient matrix $\mathbf{C}$ into the basis of molecular spinorbitals (MS), as this matrix is essential for later transformations of the core Hamiltonian and the Coulomb integrals. Let the number of spatial basis functions be denoted as $n$. The coefficient matrix in the MS basis will then be of dimension $2n \times 2n$, to incorporate the electron spin. In this transformation, each column (orbital) is duplicated and the result is vertically concatenated with itself. The corresponding segments are then zeroed out to respect the Pauli exclusion principle. This is done mathematically using the tiling matrix $\mathbf{P}_{n \times 2n}$, defined as

\begin{equation}
\mathbf{P}=\begin{pmatrix}e_1&e_1&e_2&e_2&\dots&e_n&e_n\end{pmatrix}
\end{equation}

where $e_i$ represents the $i$-th column of the identity matrix $\mathbf{I}\_n$. To facilitate the correct allocation of spin states, we define the matrices $\mathbf{M}\_{n \times 2n}$ and $\mathbf{N}_{n \times 2n}$ with elements given by

\begin{equation}
M_{ij}=1-j\bmod 2,N_{ij}=j \bmod 2
\end{equation}

The coefficient matrix $\mathbf{C}$ in the MS basis is then expressed as

\begin{equation}
\mathbf{C}^{MS}=\begin{pmatrix}\mathbf{CP}\\\ \mathbf{CP}\end{pmatrix}\odot\begin{pmatrix}\mathbf{M}\\\ \mathbf{N}\end{pmatrix}
\end{equation}

where $\odot$ denotes the Hadamard product. This transformed matrix $\mathbf{C}^{MS}$ is then used to transform the core Hamiltonian and the Coulomb integrals to the MS basis. For the core Hamiltonian, we specifically arrange the elements to be block-diagonal, with identical diagonal blocks corresponding to the spatial orbitals and zero off-diagonal blocks, mathematically written as

\begin{equation}
H_{pq}^{core,MS}=C_{\mu p}^{MS}(\mathbf{I}\_{2}\otimes_K\mathbf{H}^{core})\_{\mu\nu}C_{\nu r}^{MS}
\end{equation}

where $\otimes_K$ stands for the Kronecker product. Similarly, the transformation of the Coulomb integrals in the MS basis is given by

\begin{equation}
J_{pqrs}^{MS}=C_{\mu p}C_{\nu q}(\mathbf{I}\_{2}\otimes_K(\mathbf{I}\_{2}\otimes_K\mathbf{J})^{(4,3,2,1)})\_{\mu\nu\kappa\lambda}C_{\kappa r}C_{\lambda s}
\end{equation}

where the superscript $(4,3,2,1)$ denotes the axes transposition. This notation accounts for the spin modifications and ensures that the transformations adhere to quantum mechanical principles.

### The Determinants Generation

To proceed with the calculation, the next task involves generating all possible Slater determinants. The number of these determinants $N_D$ can be calculated using the binomial coefficients

\begin{equation}
N_D=\binom{n}{\frac{k}{2}}^2
\end{equation}

assuming $k$ is the total number of electrons, and $n$ is the total number of spinorbitals divided evenly between alpha and beta spins. Each determinant is formed by assigning the $\alpha$ and $\beta$ electrons to their respective orbitals. For practical representation, it's useful to describe determinants as arrays of numbers, where each number corresponds to the indices of the occupied orbitals. Defining the indices for potential $\alpha$ and $\beta$ orbitals, we can let $A=\\{1,3,5,\dots,n-1\\}$ represent all possible $\alpha$ orbitals and $B=\\{2,4,6,\dots,n\\}$ all possible $\beta$ orbitals. The full set of determinants, $D$, can be expressed as:

\begin{equation}
D=\\{\alpha\cup\beta\,|\,\alpha\subseteq A,\beta\subseteq B,|\alpha|=|\beta|=\frac{k}{2}\\}
\end{equation}

Here, each determinant is a combination where $\alpha$ and $\beta$ are subsets of $A$ and $B$ respectively, with each subset containing exactly $\frac{k}{2}$ orbitals, reflecting the equal distribution of electrons among $\alpha$ and $\beta$ spins.

### Slater-Condon Rules and the Eigenvalue Problem

Having generated all possible determinants, we can now proceed to compute the matrix elements of the Full Configuration Interaction (FCI) Hamiltonian, $\mathbf{H}^{CI}$. Each row and column of this matrix corresponds to a particular Slater determinant. The matrix elements between these determinants are calculated according to the Slater-Condon rules. These rules help us determine the interaction contributions by comparing two determinants and noting the number of differing spinorbitals. The computation of matrix elements $\mathbf{H}_{ij}^{CI}$ is based on how the determinants differ, and can be described by the following conditions.

\begin{equation}
\mathbf{H}_{ij}^{CI}=
\begin{cases} 
\displaystyle \sum_kH\_{kk}^{core,MS}+\frac{1}{2}\sum_k\sum_lJ\_{kkll}^{MS}-J\_{kllk}^{MS}&D_i=D_j \\\\\
\displaystyle H\_{pr}^{core,MS}+\sum_kJ\_{prkk}^{MS}-J\_{pkkr}^{MS}&D_i=\\{\dotsi p\dotsi\\}\land D_j=\\{\dotsi r\dotsi\\} \\\\\
\displaystyle \vphantom{\sum_k}J\_{prqs}^{MS}-J\_{psqr}^{MS}&D_i=\\{\dotsi p\dotsi q\dotsi\\}\land D_j=\\{\dotsi r\dotsi s\dotsi\\} \\\\\
\displaystyle \vphantom{\sum_k}0&\text{otherwise}
\end{cases}
\end{equation}

The summations extend over all orbitals common between the two compared determinants. This detailed approach ensures all electron-electron interactions are accurately accounted for in the Hamiltonian matrix, crucial for correct FCI calculations.

{:.warning}
> To accurately calculate the differences between two Slater determinants, it's essential to sometimes permute the indices of one determinant for proper alignment before comparison. For instance, consider two determinants, $D_1=\\{1,2,3,4\\}$ and $D_2=\\{1,3,2,4\\}$. To properly compare $D_2$ with $D_1$, you need to permute the indices in $D_2$ to match $D_1$, resulting in $D_2$ being reorganized as $\\{1,2,3,4\\}$. In this case, the number of differences after proper alignment is 0.
>
> Moreover, it’s crucial to remember that each permutation of the orbitals alters the sign of the determinant. Therefore, when computing the matrix elements, this change must be accounted for by multiplying the result by $(−1)^p$, where $p$ is the number of permutations made to align the determinants. This factor is vital for ensuring the correct sign of the matrix elements, which can significantly impact the results of quantum mechanical calculations.
>
> This requirement for alignment and sign adjustment adds complexity to the implementation, as each comparison between determinants not only involves checking for differences but also finding the minimal number of permutations needed for alignment.

Once the matrix elements are computed, the FCI Hamiltonian matrix $\mathbf{H}^{CI}$ is diagonalized to obtain the eigenvalues and eigenvectors. The eigenvalues represent the total energy of the system in ground and excited states, while the eigenvectors provide the coefficients for the linear combination of Slater determinants that best approximate the true wave function.

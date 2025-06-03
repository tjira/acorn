---
title: Generalized Eigenvalue Problem
parent: Mathematical Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Generalized Eigenvalue Problem<!--\label{sec:generalized_eigenvalue_problem}-->

The generalized eigenvalue problem is a mathematical problem that arises in many areas of science and engineering. It is a generalization of the standard eigenvalue problem, where we seek the eigenvalues and eigenvectors of a square matrix. In the generalized eigenvalue problem, we consider two square matrices $$\symbf{A}$$ and $$\symbf{B}$$, and we seek the eigenvalues $$\lambda$$ and eigenvectors $$\symbf{x}$$ that satisfy the equation

$$
\begin{equation}\label{eq:gen_eig}
\symbf{A}\symbf{C}=\symbf{B}\symbf{C}\symbf{\Lambda},
\end{equation}
$$

where $$\symbf{C}$$ is a matrix of eigenvectors and $$\symbf{\Lambda}$$ is a diagonal matrix of eigenvalues. The quick way to solve the generalized eigenvalue problem is to transform it into a standard eigenvalue problem by multiplying both sides of the equation by the inverse of $$\symbf{B}$$ as

$$
\begin{equation}
\symbf{B}^{-1}\symbf{A}\symbf{C}=\symbf{C}\symbf{\Lambda}.
\end{equation}
$$

This method is not always numerically stable, especially when the matrices $$\symbf{A}$$ and $$\symbf{B}$$ are ill-conditioned. Should you try to use this method for the Roothaan equations in the Hartree--Fock method, you would find that the solution is is not correct. A more stable approach is to modify the equation \eqref{eq:gen_eig} as

$$
\begin{align}
\symbf{A}\symbf{C}&=\symbf{B}\symbf{C}\symbf{\Lambda}\nonumber \\
\symbf{B}^{-\frac{1}{2}}\symbf{A}\symbf{C}&=\symbf{B}^{\frac{1}{2}}\symbf{C}\symbf{\Lambda}\nonumber \\
\symbf{B}^{-\frac{1}{2}}\symbf{A}\symbf{B}^{-\frac{1}{2}}\symbf{B}^{\frac{1}{2}}\symbf{C}&=\symbf{B}^{\frac{1}{2}}\symbf{C}\symbf{\Lambda},
\end{align}
$$

where we solve the standard eigenvalue problem for the matrix $$\symbf{B}^{-\frac{1}{2}}\symbf{A}\symbf{B}^{-\frac{1}{2}}$$ and obtain the eigenvectors $$\symbf{B}^{\frac{1}{2}}\symbf{C}$$ and eigenvalues $$\symbf{\Lambda}$$.<!--\supercite{10.48550/arXiv.1903.11240}--> The eigenvectors of the original problem are then given by $$\symbf{C}=\symbf{B}^{-\frac{1}{2}}\symbf{B}^{\frac{1}{2}}\symbf{C}$$ and the eigenvalues are the same.

{:.cite}

---
title: Generalized Eigenvalue Problem
parent: Mathematical Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Generalized Eigenvalue Problem

The generalized eigenvalue problem is a mathematical problem that arises in many areas of science and engineering. It is a generalization of the standard eigenvalue problem, where we seek the eigenvalues and eigenvectors of a square matrix. In the generalized eigenvalue problem, we consider two square matrices $\mathbf{A}$ and $\mathbf{B}$, and we seek the eigenvalues $\lambda$ and eigenvectors $\mathbf{x}$ that satisfy the equation

\begin{equation}\label{eq:gen_eig}
\mathbf{A}\mathbf{C}=\mathbf{B}\mathbf{C}\mathbf{\Lambda},
\end{equation}

where $\mathbf{C}$ is a matrix of eigenvectors and $\mathbf{\Lambda}$ is a diagonal matrix of eigenvalues. The quick way to solve the generalized eigenvalue problem is to transform it into a standard eigenvalue problem by multiplying both sides of the equation by the inverse of $\mathbf{B}$ as

\begin{equation}
\mathbf{B}^{-1}\mathbf{A}\mathbf{C}=\mathbf{C}\mathbf{\Lambda}.
\end{equation}

This method is not always numerically stable, especially when the matrices $\mathbf{A}$ and $\mathbf{B}$ are ill-conditioned. Should you try to use this method for the Roothaan equations in the Hartree--Fock method, you would find that the solution is is not correct. A more stable approach is to modify the equation \eqref{eq:gen_eig} as

\begin{align}
\mathbf{A}\mathbf{C}&=\mathbf{B}\mathbf{C}\mathbf{\Lambda}\nonumber \\\\\
\mathbf{B}^{-\frac{1}{2}}\mathbf{A}\mathbf{C}&=\mathbf{B}^{\frac{1}{2}}\mathbf{C}\mathbf{\Lambda}\nonumber \\\\\
\mathbf{B}^{-\frac{1}{2}}\mathbf{A}\mathbf{B}^{-\frac{1}{2}}\mathbf{B}^{\frac{1}{2}}\mathbf{C}&=\mathbf{B}^{\frac{1}{2}}\mathbf{C}\mathbf{\Lambda},
\end{align}

where we solve the standard eigenvalue problem for the matrix $\mathbf{B}^{-\frac{1}{2}}\mathbf{A}\mathbf{B}^{-\frac{1}{2}}$ and obtain the eigenvectors $\mathbf{B}^{\frac{1}{2}}\mathbf{C}$ and eigenvalues $\mathbf{\Lambda}$.<!--\supercite{10.48550/arXiv.1903.11240}--> The eigenvectors of the original problem are then given by $\mathbf{C}=\mathbf{B}^{-\frac{1}{2}}\mathbf{B}^{\frac{1}{2}}\mathbf{C}$ and the eigenvalues are the same.

{:.cite}

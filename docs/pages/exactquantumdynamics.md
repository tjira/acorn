---
title: Exact Quantum Dynamics
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Exact Quantum Dynamics

# Adiabatic Dynamics

\begin{equation}
\exp\left(-iH\Delta t\right)=\exp\left(-iV\frac{\Delta t}{2}\right)\exp\left(-iT\Delta t\right)\exp\left(-iV\frac{\Delta t}{2}\right)
\end{equation}

# Two-State Nonadiabatic Dynamics

\begin{equation}
\exp\left(-iH\Delta t\right)=\exp\left(-i\begin{pmatrix}V_{11}&V_{12}\\\V_{21}&V_{22}\end{pmatrix}\frac{\Delta t}{2}\right)\exp\left(-i\begin{pmatrix}T&0\\\0&T\end{pmatrix}\Delta t\right)\exp\left(-i\begin{pmatrix}V_{11}&V_{12}\\\V_{21}&V_{22}\end{pmatrix}\frac{\Delta t}{2}\right)
\end{equation}

\begin{equation}
\exp\left[-i\begin{pmatrix}V_{11}&V_{12}\\\V_{21}&V_{22}\end{pmatrix}\frac{\Delta t}{2}\right]=\exp\left[-i(V_{11}+V_{22})\frac{\Delta t}{4}\right]\left[\cos\left(\sqrt{D}\frac{\Delta t}{4}\right)\begin{pmatrix}1&0\\\0&1\end{pmatrix}+i\frac{\sin\left(\sqrt{D}\frac{\Delta t}{4}\right)}{\sqrt{D}}\begin{pmatrix}V_{22}-V_{11}&-2V_{12}\\\ -2V_{21}&V_{11}-V_{22}\end{pmatrix}\right],
\end{equation}
where $D=4|V_{21}|^2+(V_{11}-V_{22})^2$.


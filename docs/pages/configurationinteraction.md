---
title: Configuration Interaction
parent: Electronic Structure Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Configuration Interaction

Case 1: $\ket{K}=\ket{\dotsi mn\dotsi}$
\begin{equation}
\braket{K|H|K}=\sum_m^N[m|h|m]+\frac{1}{2}\sum_m^N\sum_n^N[mm|nn]-[mn|nm]
\end{equation}

Case 2: $\ket{K}=\ket{\dotsi mn\dotsi}$, $\ket{L}=\ket{\dotsi pn\dotsi}$
\begin{equation}
\braket{K|H|L}=[m|h|p]+\sum_m^N[mp|nn]-[mn|np]
\end{equation}

Case 3: $\ket{K}=\ket{\dotsi mn\dotsi}$, $\ket{L}=\ket{\dotsi pq\dotsi}$
\begin{equation}
\braket{K|H|L}=[mp|nq]-[mq|np]
\end{equation}

The brakets in the table above are all in chemists' notation.

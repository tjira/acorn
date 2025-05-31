---
title: Electronic Structure Methods
has_children: true
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Electronic Structure Methods<!--\label{sec:electronic_structure_methods}-->

This part provides an educational exploration into the computational techniques fundamental to understanding molecular electronic structure in quantum chemistry. Beginning with the Hartree--Fock method, the text introduces this foundational approach for determining molecular orbitals and electronic energies by approximating the interactions of electrons through a mean-field approximation. The HF method forms the basis for subsequent methods and is presented with a practical coding exercise that guides readers in implementing and calculating Hartree--Fock energies in Python.

Moving beyond HF, the text delves into Møller--Plesset Perturbation Theory, which improves Hartree--Fock predictions by introducing corrections for electron correlation through a perturbative approach. This section includes exercises on calculating second- and third-order corrections, allowing readers to enhance their understanding of how electron interactions can be more accurately incorporated. Configuration Interaction theory is then presented as an approach for representing the molecular wavefunction as a combination of electron configurations. Here, readers learn the theoretical basis of Configuration Interaction and engage with practical examples focused on constructing the Configuration Interaction Hamiltonian matrix and solving for molecular energies, particularly emphasizing Full Configuration Interaction for high accuracy in small systems.

The part culminates with a discussion on the Coupled Cluster theory, a highly accurate and computationally efficient method for capturing electron correlation effects, often used for small to medium-sized systems. By introducing truncations such as Coupled Cluster Doubles and Coupled Cluster Singles and Doubles, the text demonstrates how electron correlation can be systematically included while balancing computational cost. The Coupled Cluster section provides iterative coding exercises for calculating correlation energies, rounding out the document’s comprehensive approach to electronic structure methods in computational quantum chemistry. Through this blend of theory, mathematical formulations, and hands-on coding exercises, the document serves as an invaluable resource for building a strong foundational understanding of electronic structure methods.

{:.no_toc .text-delta}

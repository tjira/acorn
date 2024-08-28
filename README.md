<h1 align="center">Quantum Acorn</h1>

<h4 align="center">
  <a href="https://github.com/tjira/acorn#features">Features</a>
  ·
  <a href="https://github.com/tjira/acorn#compilation">Compilation</a>
  ·
  <a href="https://github.com/tjira/acorn#examples">Examples</a>
  ·
  <a href="https://github.com/tjira/acorn#credits">Credits</a>
  ·
  <a href="https://tjira.github.io/acorn/">Docs</a>
</h4>

<p align="center">
    <a href="https://github.com/tjira/acorn/pulse">
        <img src="https://img.shields.io/github/last-commit/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/blob/master/LICENSE.md">
        <img src="https://img.shields.io/github/license/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/actions/workflows/test.yml">
        <img src="https://img.shields.io/github/actions/workflow/status/tjira/acorn/test.yml?style=for-the-badge&label=test"/>
    </a>
    <a href="https://app.codecov.io/gh/tjira/acorn">
        <img src="https://img.shields.io/codecov/c/github/tjira/acorn?style=for-the-badge"/>
    </a>
    <br>
    <a href="https://github.com/tjira/acorn/stargazers">
        <img src="https://img.shields.io/github/stars/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn">
        <img src="https://img.shields.io/github/languages/code-size/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/releases/latest">
        <img src="https://img.shields.io/github/v/release/tjira/acorn?display_name=tag&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/releases/latest">
        <img src="https://img.shields.io/github/downloads/tjira/acorn/total?style=for-the-badge"/>
    </a>
</p>

<p align="center">
Quantum Acorn, a collection of electronic structure methods compiled into a dependency-free binaries. If you are here for the educational scripts, you can find them in the education folder.
</p>

## Features

Below are all the important features of Acorn divided into categories. If you are looking for the educational scripts, you can find them in the education folder.

### Quantum Mechanical Methods

* Numerically Exact Adiabatic & Nonadiabatic Quantum Dynamics
* Hartree-Fock Method & Møller–Plesset Perturbation Theory
* Configuration Interaction & Coupled Clusters Methods

## Compilation

Before the compilation process, make sure you have [eigen](https://gitlab.com/libeigen/eigen), [boost](https://github.com/boostorg/boost) and [openblas](https://github.com/OpenMathLib/OpenBLAS) libraries installed. On debian-based distributions, you can do it with the following command.

```bash
sudo apt install libboost-dev libeigen3-dev libopenblas-dev
```

After these well-estabilished libraries are installed, you need to compile several other dependencies. Fortunately, the Acorn package contains a script to install them all-in-one. You can install all additional dependencies by executing `./script/liball.sh` from the root directory. Now, go grab some coffee, it can take several hours. After the library compilation is finished, the project configuration should finish without errors.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
```

And we can build with the following command.

```bash
cmake --build build
```

After the compilation the bin folder will be created along with the executables.

## Examples

All the examples are implemented as `example/makefile` targets. If you are in example directory, you can run the HF example on a water molecule simply run the `make hf` command. The calculation should finish without errors. Feel free to explore all the examples.

## Credits

* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
* [eigen](https://gitlab.com/libeigen/eigen) - C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [exprtk](https://github.com/ArashPartow/exprtk) - C++ Mathematical Expression Parsing and Evaluation Library.
* [fftw](https://www.fftw.org) - C Subroutine Library for Computing the Discrete Fourier Transform.
* [libint](https://github.com/evaleev/libint) - High-Performance Library for Computing Gaussian Integrals in Quantum Mechanics.
* [pytorch](https://github.com/pytorch/pytorch) - Tensors and Dynamic neural networks in Python with strong GPU acceleration.

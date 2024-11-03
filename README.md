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
Quantum Acorn, a collection of electronic structure methods compiled into a dependency-free binary. If you are here for the educational scripts, you can find them in the education folder.
</p>

## Features

Below are all the important features of Acorn divided into categories.

### Quantum Mechanical Methods

* Numerically Exact Adiabatic & Nonadiabatic Quantum Dynamics with Arbitrary Number of States & Dimensions
* Restricted Hartree–Fock for Closed Shell & Generalized Hartree–Fock for Open Shell Systems
* Møller–Plesset Perturbation Theory, Configuration Interaction & Coupled Cluster Methods

## Compilation

All the libraries the program needs are included in the installation script. To install all the libraries, navigate to the project root directory and execute the following command.

```bash
./script/general/library.sh SHARED 1
```

You can also perform a static compilation. The number at the end is number of cores the compilation process uses. The program needs some heavy libraries so it takes some time. After the library compilation finishes, all the header files and compiled libraries can be found in the *external* directory. If you don't want to wait for the compilation process, you can also download libraries used in the latest release by running the following command.

```bash
./script/general/libdown.sh SHARED
```

This command also creates the *external* directory with libraries. After you have obtained the compiled libraries, configure the project by the following command.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DSTATIC=OFF
```

And build with the following command.

```bash
cmake --build build
```

After the compilation, the bin folder will be created along with the executable.

## Examples

All the examples are located in `example/input` folder. To run an example Hartree-Fock calculation, execute the corresponding example file using the `acorn -i example/input/hf.json` command from the project root. Feel free to explore all the examples.

## Credits

* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
* [eigen](https://gitlab.com/libeigen/eigen) - C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [exprtk](https://github.com/ArashPartow/exprtk) - C++ Mathematical Expression Parsing and Evaluation Library.
* [fftw](https://www.fftw.org) - C Subroutine Library for Computing the Discrete Fourier Transform.
* [json](https://github.com/nlohmann/json) - JSON for Modern C++.
* [libint](https://github.com/evaleev/libint) - High-Performance Library for Computing Gaussian Integrals in Quantum Mechanics.
* [pytorch](https://github.com/pytorch/pytorch) - Tensors and Dynamic neural networks in Python with strong GPU acceleration.

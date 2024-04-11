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
    <a href="https://github.com/tjira/acorn/stargazers">
        <img src="https://img.shields.io/github/stars/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/releases/latest">
        <img src="https://img.shields.io/github/downloads/tjira/acorn/total?style=for-the-badge"/>
    </a>
    <br>
    <a href="https://github.com/tjira/acorn">
        <img src="https://tokei.rs/b1/github/tjira/acorn?category=code&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn">
        <img src="https://img.shields.io/github/languages/code-size/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/releases/latest">
        <img src="https://img.shields.io/github/v/release/tjira/acorn?display_name=tag&style=for-the-badge"/>
    </a>
</p>

<p align="center">
Quantum Acorn, a collection of electronic structure methods compiled into a dependency-free binary. If you are here for the educational scripts, you can find them in the education folder.
</p>

## Features

Below are all the important features of Acorn divided into categories. If you are looking for the educational scripts, you can find them in the education folder.

### Quantum Mechanical Methods

* Numerically Exact 1D and 2D Adiabatic Quantum Dynamics
* Numerically Exact 1D Nonadiabatic Quantum Dynamics
* Hartree-Fock Method (RHF & UHF)
* Møller–Plesset Perturbation Theory
* Configuration Interaction (FCI, CISD, CID, CIS)

### Additional Calculations

* Analytical Gradient for RHF and Numerical for Everything Else
* Numerical Hessians and Frequency Analysis for Every Method
* Molecular Dynamics With Calculated or Provided Gradients
* Mulliken Population Analysis for RHF
* Berendsen Thermostatting for MD

## Compilation

The software requires the [libint](https://github.com/evaleev/libint) library. Before the compilation process, make sure you have [eigen](https://gitlab.com/libeigen/eigen) and [boost](https://github.com/boostorg/boost) installed. On debian-based distributions, you can do it with the following command.

```bash
sudo apt install libboost-dev libeigen3-dev
```

To compile the library execute `./script/libint.sh` from the project root directory. This command creates the `libint` directory with the compiled library. Now, we export the necessary environment variables.

```bash
export CPLUS_INCLUDE_PATH="$PWD/libint/install/include:$CPLUS_INCLUDE_PATH"
export LIBRARY_PATH="$PWD/libint/install/lib:$LIBRARY_PATH"
```

After this, the project configuration should finish without errors.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
```

And we can build with the following command.

```bash
cmake --build build
```

After the compilation the bin folder will be created along with the executable.

## Examples

All the example inputs are located in the `example/input` folder. They are meant to be kept there due to the relative paths to the molecules. If you are in the project root directory, you can run one of the examples with the following command.

```bash
./bin/acorn example/input/rhf.json
```

The calculation should finish without errors. Feel free to explore all the examples. Keep in mind that to execute the ORCA dynamics example you need the ORCA executable in your PATH variable.

## Credits

* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
* [exprtk](https://github.com/ArashPartow/exprtk) - C++ Mathematical Expression Parsing and Evaluation Library.
* [fftw](https://www.fftw.org) - C Subroutine Library for Computing the Discrete Fourier Transform .
* [json](https://github.com/nlohmann/json) - JSON for Modern C++.
* [libint](https://github.com/evaleev/libint) - High-Performance Library for Computing Gaussian Integrals in Quantum Mechanics.

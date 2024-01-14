<h1 align="center">Quantum Acorn</h1>

<h4 align="center">
  <a href="https://github.com/tjira/acorn#compilation">Compilation</a>
  ·
  <a href="https://github.com/tjira/acorn#examples">Examples</a>
  ·
  <a href="https://github.com/tjira/acorn#printing">Printing</a>
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
        <img src="https://img.shields.io/github/v/release/tjira/acorn?display_name=tag&style=for-the-badge"/>
    </a>
    <br>
    <a href="https://github.com/tjira/acorn">
        <img src="https://tokei.rs/b1/github/tjira/acorn?category=code&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn">
        <img src="https://img.shields.io/github/languages/code-size/tjira/acorn?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/acorn/releases/latest">
        <img src="https://img.shields.io/github/downloads/tjira/acorn/total?style=for-the-badge"/>
    </a>
</p>

<p align="center">
Quantum Acorn, a dynamic collection of electronic structure methods, effortlessly transforms input geometry into quantum insights. Simply pass an input file, and watch as results appear in your terminal.
</p>

## Features

Below are all the important features of Acorn divided into categories.

### Quantum Mechanical Methods

* Hartree-Fock Method
* Møller–Plesset Perturbation Theory

### Additional Calculations

* Gradients, Hessians and Frequency Analysis

## Compilation

The software requires the [libint](https://github.com/evaleev/libint) library. Before the library compilation process, make sure you have [eigen](https://gitlab.com/libeigen/eigen) and [boost](https://github.com/boostorg/boost) installed. On debian-based distributions, you can do it with the following command.

```bash
sudo apt install libboost-all-dev libeigen3-dev
```

To compile the library execute `./script/libint.sh` from the project root directory. This command creates the *libint* folder with the compiled library. Now, we export the necessary environment variables.

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

## Credits

* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
* [json](https://github.com/nlohmann/json) - JSON for Modern C++.
* [libint](https://github.com/evaleev/libint) - High-performance library for computing Gaussian integrals in quantum mechanics.

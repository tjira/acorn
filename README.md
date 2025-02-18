<h1 align="center">Quantum Acorn</h1>

<h4 align="center">
  <a href="https://github.com/tjira/acorn#features">Features</a>
  ·
  <a href="https://github.com/tjira/acorn#compilation">Compilation</a>
  ·
  <a href="https://github.com/tjira/acorn#examples">Examples</a>
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
    <a href="https://github.com/tjira/acorn/stargazers">
        <img src="https://img.shields.io/github/stars/tjira/acorn?style=for-the-badge"/>
    </a>
    <br>
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
* Restricted Hartree–Fock Method, Møller–Plesset Perturbation Theory & Configuration Interaction Methods

## Compilation

Since the software is coded in zig, you need to have the zig compiler installed on your system. You can download the latest version of the zig compiler from the [official website](https://ziglang.org/download). After you have installed the zig compiler, navigate to the project root and run the following command to compile the project.

```bash
zig build --release=fast --summary all
```

This will compile the project and create a binary file named `acorn` in the `zig-out/arch-os` folder, where `arch` is the architecture of your system and `os` is the operating system you are using. If you don't know your CPU architecture, you probably want the `x86_64` binary. You can also perform tests on the project by running the following command.

```bash
zig build --release=fast --summary all test
```

If some tests fail, let me know by creating an issue. If all the tests pass, you can run the binary file using the following command.

```bash
./zig-out/arch-os/acorn
```

You should see the version of the compiler and execution time of the program. If you see this, the program is working correctly.

## Examples

Below are some examples of the quantum mechanical methods implemented in Acorn. The plotting commands use the script `plot.py` in the `python` folder. You can use this script to plot the data generated by the program.


### Real Time Quantum Dynamics

This example demonstrates the real-time quantum dynamics of a particle in a harmonic potential. The particle is initialized in a Gaussian wavepacket and allowed to evolve in time. The wavefunction is written to a file named `WAVEFUNCTION.mat` in the current directory. The input file for this example is shown below.

```json
{
    "quantum_dynamics" : {
        "adiabatic" : false,
        "iterations" : 500,
        "mode" : [0, 1],
        "time_step" : 0.1,
        "grid" : {
            "limits" : [-8, 8],
            "points" : 1024
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0],
            "position" : [1],
            "state" : 0,
            "gamma" : 2
        },
        "write" : {
            "wavefunction" : "WAVEFUNCTION.mat"
        },
        "potential" : "harmonic1D_1"
    }
}
```

The input file can be run like any other program in Acorn, no special flags are required. This simulation is fast and should complete in under a second. After the simulation is complete, you can visualize the results using the commands in the below sections.

<details> <summary><b>Wavefunction</b></summary>

```bash
plot.py WAVEFUNCTION.mat:0,1 --animate 2 --xlabel "Coordinate (a.u.)" --ylabel "Wavefunction"
```

![Wavefunction](graphics/rtpa_wavefunction.gif)
</details>

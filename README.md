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
Quantum Acorn, a collection of electronic structure methods compiled into a small binary. If you are here for the educational scripts, you can find them in the education folder.
</p>

## Features

Below are all the important features of Acorn divided into categories.

### Time-Independent Quantum Mechanics

* Restricted Hartree–Fock Method with DIIS
* Møller–Plesset Perturbation Theory
* Full Configuration Interaction

### Time-Dependent Quantum Mechanics

* Quantum Dynamics with Arbitrary Number of States & Dimensions
* Fewest Switches Surface Hopping Dynamics
* Landau–Zener Surface Hopping Dynamics
* Mapping Approach to Surface Hopping Dynamics

## Compilation

Before the compilation, make sure you have `cblas`, `fftw`, `gsl` and `lapacke` libraries installed on your system. You can install them using your package manager. For example, on Ubuntu, you can run the following command.

```bash
sudo apt install fftw3 libatlas-base-dev libgsl-dev liblapacke-dev
```

After the necessary dependencies you need to have the zig compiler installed on your system. You can download the latest version of the zig compiler from the [official website](https://ziglang.org/download) or install it from your favourite package manager. After you have installed the zig compiler, navigate to the project root and run the following command to compile the project.

```bash
zig build && zig build script
```

This will compile the project and create a binary file named `acorn` in the `zig-out/arch-os` folder, where `arch` is the architecture of your system and `os` is the operating system you are using. It will also generate some usefull wrappers around the binary in the same folder. You can also perform tests on the project by running the following command.

```bash
zig build test
```

If some tests fail, let me know by creating an issue. If all the tests pass, you can run the binary file using the following command.

```bash
./zig-out/arch-os/acorn
```

You should see the version of the compiler and execution time of the program. If you see this, the program is working correctly.

## Examples

Below are some examples of the quantum mechanical methods implemented in Acorn. The plotting commands use the script `plot.py` in the `python` folder. You can use this script to plot the data generated by the program.

- [Hartree–Fock Method](#hartreefock-method)
- [Møller–Plesset Perturbation Theory](#møllerplesset-perturbation-theory)
- [Real Time Adiabatic Quantum Dynamics](#real-time-adiabatic-quantum-dynamics)
- [Real Time Nonadiabatic Quantum Dynamics](#real-time-nonadiabatic-quantum-dynamics)
- [Surface Hopping Dynamics](#surface-hopping-dynamics)

### Hartree–Fock Method

The following input calculates the Hartree–Fock energy of a provided system. The system in the example is somewhat optimized water. If you want to provide an .xyz file with the coordinates of the system, simply delete the `system` field and add the `system_file` field with the path to the geometry file.

```json
{
    "hartree_fock" : {
        "system" : {
            "atoms" : [8, 1, 1],
            "coords" : [
                [-0.04, -0.01, -0.01],
                [ 0.65, -0.51,  0.47],
                [ 0.38,  0.87, -0.06]
            ]
        },
        "integral" : {
            "basis" : "sto-3g"
        }
    }
}
```

### Møller–Plesset Perturbation Theory

The following input calculates energy of a system using the MP2 method. The system is the same as for the Hartree–fock method.

```json
{
    "moller_plesset" : {
        "hartree_fock" : {
            "system" : {
                "atoms" : [8, 1, 1],
                "coords" : [
                    [-0.04, -0.01, -0.01],
                    [ 0.65, -0.51,  0.47],
                    [ 0.38,  0.87, -0.06]
                ]
            },
            "integral" : {
                "basis" : "sto-3g"
            }
        },
        "order" : 2
    }
}
```

### Real Time Adiabatic Quantum Dynamics

This example demonstrates the real-time quantum dynamics of a particle in a harmonic potential. The particle is initialized in a Gaussian wavepacket and allowed to evolve in time. The wavefunction is written to a file named `WAVEFUNCTION.mat` in the current directory. The input file for this example is shown below.

```json
{
    "quantum_dynamics" : {
        "adiabatic" : false,
        "iterations" : 1000,
        "mode" : [0, 1],
        "time_step" : 0.1,
        "grid" : {
            "limits" : [-8, 8],
            "points" : 512
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0],
            "position" : [1],
            "state" : 0,
            "gamma" : 2
        },
        "write" : {
            "wavefunction" : "WAVEFUNCTION.mat",
            "autocorrelation_function" : "ACF.mat",
            "spectrum" : "SPECTRUM.mat"
        },
        "potential" : "harmonic1D_1"
    }
}
```

The input file can be run like any other program in Acorn, no special flags are required. This simulation is fast and should complete in under a second. You can visualize the wavefunction, autocorrelation function or the vibrational spectrum with the commands below.

```bash
lines.py WAVEFUNCTION.mat:0,1 --legends every "Re(\$\Psi_0\$)" "Im(\$\Psi_0\$)" --xlabel "Coordinate (a.u.)" --ylabel "Wavefunction" --animate 2
```

```bash
lines.py ACF.mat:0,1 SPECTRUM.mat --figsize 6 16 --legends 0,1 "Re(\$<\Psi_0|\Psi>\$)" "Im(\$<\Psi_0|\Psi>\$)" --subplots 121 121 122 --xlabel "Time (a.u.)" "Energy (a.u.)" --xlim nan nan 0 6 --ylabel "ACF" "Intensity"
```

Acorn also supports higher dimensions. As an example you can simulate a 2D wavefunction in a 2D harmonic potential using the following input.

```json
{
    "quantum_dynamics" : {
        "adiabatic" : false,
        "iterations" : 1000,
        "mode" : [0, 1],
        "time_step" : 0.1,
        "grid" : {
            "limits" : [-8, 8],
            "points" : 256
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0, 0],
            "position" : [1, 1],
            "state" : 0,
            "gamma" : 2
        },
        "write" : {
            "wavefunction" : "WAVEFUNCTION.mat",
            "autocorrelation_function" : "ACF.mat",
            "spectrum" : "SPECTRUM.mat"
        },
        "potential" : "harmonic2D_1"
    }
}
```

This simulation takes a few seconds, since the time complexity increases exponentially. Visualizing the 3D complex wavefunction is a little tricky. One way is to plot the square of the wavefunction on a 2D heatmap. You can visualize the wavefunction this way and ACF with spectrum the same way as above with the commands below.

```bash
heatmap.py WAVEFUNCTION.mat:0,1 --xlabel "Coordinate (a.u.)" --ylabel "Coordinate (a.u.)" --transform norm --animate 2
```

```bash
lines.py ACF.mat:0,1 SPECTRUM.mat --figsize 6 16 --legends 0,1 "Re(\$<\Psi_0|\Psi>\$)" "Im(\$<\Psi_0|\Psi>\$)" --subplots 121 121 122 --xlabel "Time (a.u.)" "Energy (a.u.)" --xlim nan nan 0 9 --ylabel "ACF" "Intensity"
```

### Real Time Nonadiabatic Quantum Dynamics

This example demonstrates the real-time quantum dynamics of a first Tully potential in the adiabatic basis. It exports the resulting time-dependent wavefunction and the potential. The input file for this example is shown below.

```json
{
    "quantum_dynamics" : {
        "adiabatic" : true,
        "iterations" : 500,
        "mode" : [0, 1],
        "time_step" : 10,
        "grid" : {
            "limits" : [-16, 32],
            "points" : 2048
        },
        "initial_conditions" : {
            "mass" : 2000,
            "momentum" : [10],
            "position" : [-10],
            "state" : 1,
            "gamma" : 2
        },
        "write" : {
            "wavefunction" : "WAVEFUNCTION.mat",
            "population" : "POPULATION.mat"
        },
        "potential" : "tully1D_1"
    }
}
```

This simulation is fast and should complete approximately in a second. You can visualize the wavefunction and state population with the commands below. The wavefunctions on each states are vertically separated for better visualization.

```bash
lines.py WAVEFUNCTION.mat:0,1,2,3 --legends every "Re(\$\Psi_0\$)" "Im(\$\Psi_0\$)" "Re(\$\Psi_1\$)" "Im(\$\Psi_1\$)" --offsets every -1 -1 1 1 --xlabel "Coordinate (a.u.)" --ylabel "Wavefunction" --animate 4
```

```bash
lines.py POPULATION.mat --legends every "S\$_0\$" "S\$_1\$" --xlabel "Time (a.u.)" --ylabel "Population"
```

### Surface Hopping Dynamics

This example demonstrates how to run a surface hopping dynamics. The below example executes a Fewest Switches Surface Hopping (FSSH) dynamics on the same potential as the real-time nonadiabatic dynamics example above.

```json
{
    "classical_dynamics" : {
        "adiabatic" : true,
        "iterations" : 5000,
        "time_step" : 1,
        "trajectories" : 1000,
        "initial_conditions" : {
            "mass" : [2000],
            "momentum_mean" : [10],
            "momentum_std" : [1],
            "position_mean" : [-10],
            "position_std" : [0.5],
            "state" : [0, 1]
        },
        "log_intervals" : {
            "iteration" : 500,
            "trajectory" : 100
        },
        "write" : {
            "population_mean" : "POPULATION_MEAN.mat"
        },
        "potential" : "tully1D_1",
        "fewest_switches" : {}
    }
}
```

This simulation is slow and will take a few second to complete. You can visualize the mean population of each state with the command below.

```bash
lines.py POPULATION_MEAN.mat --legends every "S\$_0\$" "S\$_1\$" --xlabel "Time (a.u.)" --ylabel "Population"
```

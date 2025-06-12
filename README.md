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
Quantum Acorn, a collection of electronic structure methods compiled into a small dependency-free binary. If you are here for the educational scripts, you can find them in the education folder.
</p>

## Features

Acorn offers a comprehensive suite of tools for both time-independent and time-dependent quantum mechanical simulations. Key features are grouped below for clarity.

### Time-Independent Quantum Mechanics

* **Hartree–Fock Methods**  
  Includes both restricted (RHF) and generalized (GHF) Hartree–Fock, enhanced with DIIS convergence acceleration.

* **Post-Hartree–Fock Methods**  
  Support for Møller–Plesset perturbation theory (MP2) and configuration interaction (CASCI) techniques.

* **Analytical Tools**  
  Computation of energy gradients, Hessians, vibrational frequencies, and geometry optimizations for all supported methods.

### Time-Dependent Quantum Mechanics

* **Quantum Dynamics Simulations**  
  Simulate quantum dynamics across arbitrary dimensions and electronic states using customizable potential energy surfaces.

* **Bohmian Dynamics**  
  Propagate trajectory ensembles guided by the time-dependent Schrödinger equation for phase-aware quantum dynamics.

* **Surface Hopping Algorithms**  
  Perform nonadiabatic dynamics with support for FSSH, LZSH, MASH, and κTSH methods, all on user-defined potentials.

## Acquiring the Software

### Downloading Release

You can download the latest release [here](https://github.com/tjira/acorn/releases/latest). Two versions of the executable are available, with one built using musl libc and the other using GNU libc. The musl version is a statically linked ELF executable that embeds its own C standard library and has no external runtime dependencies, allowing it to run on any Linux kernel regardless of the user space environment. The GNU version, linked against GNU libc, depends on the host system's runtime libraries but can offer better integration with glibc-based distributions.

### Compilation

To install all necessary prerequisites, you need to have a few libraries installed on your system. You can download and compile them in the project root using the practical `./script/library.sh` script. The script, among other things, downloads the zig compiler into the `zig-bin` directory. After the libraries have compiled, navigate to the project root and run the following command to compile the project.

```bash
./zig-bin/zig build && zig build script
```

This will compile the project and create a binary file named `acorn` in the `zig-out/arch-os-abi` folder, where `arch` is the architecture of your system, `os` is the operating system you are using and the default ABI is musl. It will also generate some usefull wrappers around the binary in the same folder. You can also perform tests on the project by running the following command.

```bash
./zig-bin/zig build test
```

If some tests fail, let me know by creating an issue. If all the tests pass, you can run the binary file using the following command.

```bash
./zig-out/arch-os-abi/acorn
```

You should see the version of the compiler and execution time of the program. If you see this, the program is working correctly.

## Examples

Below are some examples of the quantum mechanical methods implemented in Acorn divided in categories. The plotting commands use the script `plot.py` in the `python` folder. You can use this script to plot the data generated by the program.

<details>

<summary><b>Electronic Structure Methods</b></summary>

- [Hartree–Fock Method](#hartreefock-method)
- [Møller–Plesset Perturbation Theory](#møllerplesset-perturbation-theory)
- [Full Configuration Interaction](#full-configuration-interaction)

### Hartree–Fock Method

The input below performs a geometry optimization using the Hartree–Fock method, followed by calculations of the energy, gradient, Hessian, and related properties. The example system used here is a partially optimized water molecule. To use a custom geometry from an .xyz file, remove the `system` field and add a `system_file` field pointing to the path of your coordinate file.

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
        },
        "optimize" : {},
        "gradient" : {},
        "hessian" : {}
    }
}
```

### Møller–Plesset Perturbation Theory

The input below performs a geometry optimization using the Møller–Plesset perturbation theory (MP2) method, followed by calculations of the energy, gradient, Hessian, and related properties. The example system is the same partially optimized water molecule used in the Hartree-Fock example above.

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
        "order" : 2,
        "optimize" : {},
        "gradient" : {},
        "hessian" : {}
    }
}
```

### Full Configuration Interaction

The following example calculates the Full Configuration Interaction (FCI) energy of the same water molecule used in the previous examples. You can perform geometry optimization or compute gradients and Hessians in the same way as for the Hartree–Fock (HF) or Møller–Plesset (MP2) methods. If the `active_space` variable is set to `null`, the program uses the full FCI method. Alternatively, you can define an active space by providing an array of two numbers: the first specifies the number of electrons, and the second the number of spin orbitals they occupy.

```json
{
    "configuration_interaction" : {
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
        "active_space" : null
    }
}
```

</details>

<details>

- [Real Time Adiabatic Quantum Dynamics](#real-time-adiabatic-quantum-dynamics-with-bohmian-trajectories)
- [Real Time Nonadiabatic Quantum Dynamics](#real-time-nonadiabatic-quantum-dynamics)

<summary><b>Quantum Dynamics</b></summary>

### Real Time Adiabatic Quantum Dynamics with Bohmian Trajectories

This example demonstrates real-time quantum dynamics of a particle in a harmonic potential, along with Bohmian dynamics for quantum trajectories. The particle is initially prepared as a Gaussian wavepacket and evolves over time. The resulting wavefunction is saved to a file named `WAVEFUNCTION.mat` in the current directory. The input file for this example is shown below.

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
        "hamiltonian" : {
            "dims" : 1,
            "matrix" : [["0.5*r1^2"]]
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0],
            "position" : [1],
            "state" : 0,
            "gamma" : 2
        },
        "write" : {
            "bohm_position" : "POSITION_BOHM.mat",
            "bohm_position_mean" : "POSITION_BOHM_MEAN.mat",
            "position" : "POSITION_EXACT.mat",
            "wavefunction" : "WAVEFUNCTION.mat"
        },
        "bohmian_dynamics" : {
            "trajectories" : 100
        }
    }
}
```

The input file can be run like any other program in Acorn, no special flags are required. This simulation is fast and should complete in under a second. You can visualize the wavefunction and the trajectories with the commands below.

```bash
./python/lines.py WAVEFUNCTION.mat:0,1 --legends every "Re(\$\Psi\$)" "Im(\$\Psi\$)" --xlabel "Coordinate (a.u.)" --ylabel "Wavefunction" --animate 2
```

```bash
./python/lines.py POSITION_BOHM.mat POSITION_EXACT.mat POSITION_BOHM_MEAN.mat --alphas 0-99 0.1 --colors 0-99,100,101 tab:blue tab:orange black --styles 101 dashed --xlabel "Time (a.u.)" --ylabel "Coordinate (a.u.)"
```

Acorn also supports higher dimensions. As an example you can simulate a 2D wavefunction in a 2D harmonic potential using the following input. The Bohmian dynamics is also available for higher dimension the same way as in the previous example.

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
        "hamiltonian" : {
            "dims" : 2,
            "matrix" : [["0.5*(r1^2+r2^2)"]]
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0, 0],
            "position" : [1, 1],
            "state" : 0,
            "gamma" : 2
        },
        "write" : {
            "wavefunction" : "WAVEFUNCTION.mat"
        }
    }
}
```

This simulation takes a few seconds, since the time complexity increases exponentially. Visualizing the 3D complex wavefunction is a little tricky. One way is to plot the square of the wavefunction on a 2D heatmap. You can visualize the wavefunction this way with the command below.

```bash
./python/heatmap.py WAVEFUNCTION.mat:0,1 --xlabel "Coordinate (a.u.)" --ylabel "Coordinate (a.u.)" --transform norm --animate 2
```

### Real Time Nonadiabatic Quantum Dynamics

This example illustrates real-time quantum dynamics on the first Tully potential, represented in the adiabatic basis. While the Tully potentials are predefined, you may also supply a custom potential matrix, as shown in earlier examples. Bohmian trajectories can be enabled in the same way as before. The simulation exports the time-dependent wavefunction and the potential. The input file for this example is shown below.

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
        "hamiltonian" : {
            "name" : "tully1D_1"
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
        }
    }
}
```

This simulation is fast and should complete approximately in a second. You can visualize the wavefunction and state population with the commands below. The wavefunctions on each states are vertically separated for better visualization.

```bash
./python/lines.py WAVEFUNCTION.mat:0,1,2,3 --legends every "Re(\$\Psi_0\$)" "Im(\$\Psi_0\$)" "Re(\$\Psi_1\$)" "Im(\$\Psi_1\$)" --offsets every -1 -1 1 1 --xlabel "Coordinate (a.u.)" --ylabel "Wavefunction" --animate 4
```

```bash
./python/lines.py POPULATION.mat --legends every "S\$_0\$" "S\$_1\$" --xlabel "Time (a.u.)" --ylabel "Population"
```

</details>

<details>

<summary><b>Mixed Quantum-Classical Dynamics</b></summary>

- [Surface Hopping Dynamics](#surface-hopping-dynamics)

### Surface Hopping Dynamics

This example demonstrates how to run a surface hopping dynamics. The below example executes a Fewest Switches Surface Hopping (FSSH) dynamics on the same potential as the real-time nonadiabatic quantum dynamics example above.

```json
{
    "classical_dynamics" : {
        "iterations" : 5000,
        "time_step" : 1,
        "trajectories" : 1000,
        "hamiltonian" : {
            "name" : "tully1D_1"
        },
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
        "fewest_switches" : {}
    }
}
```

This simulation is slow and will take a few second to complete. You can visualize the mean population of each state with the command below.

```bash
./python/lines.py POPULATION_MEAN.mat --legends every "S\$_0\$" "S\$_1\$" --xlabel "Time (a.u.)" --ylabel "Population"
```

</details>

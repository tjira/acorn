//! Module where all the output structs are stored.

const std = @import("std");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The atomic integral options
pub fn AtomicIntegralOptions(comptime T: type) type {
    return struct {
        pub const Print = struct {
            overlap: bool = false,
            kinetic: bool = false,
            nuclear: bool = false,
            coulomb: bool = false,
        };
        pub const Write = struct {
            overlap: ?[]const u8 = null,
            kinetic: ?[]const u8 = null,
            nuclear: ?[]const u8 = null,
            coulomb: ?[]const u8 = null,
        };

        system_file: ?[]const u8 = null,
        basis:        []const u8,

        print: Print = .{}, system: HartreeFockOptions(T).System = .{}, write: Write = .{}
    };
}

/// The classical dynamics options
pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            position_mean:          []const T,
            position_std:           []const T,
            momentum_mean:          []const T,
            momentum_std:           []const T,
            state:                  []const T,
            mass:                   []const T
        };
        pub const FewestSwitches = struct {
            quantum_substep: u32 = 10,
            decoherence_alpha: ?T = null
        };
        pub const Hamiltonian = struct {
            dims: ?u32 = null, matrix: ?[]const []const []const u8 = null, name: ?[]const u8 = null, file: ?[]const u8 = null
        };
        pub const LandauZener = struct {
            mode: ?[]const u8 = "maxdiff"
        };
        pub const SpinMapping = struct {
            pub const Resample = struct {
                reflect: bool = false, threshold: T = 0
            };

            resample: ?Resample = null
        };
        pub const LogIntervals = struct {
            trajectory: u32 = 1, iteration: u32 = 1
        };
        pub const Write = struct {
            bloch_vector:                  ?[]const u8 = null,
            bloch_vector_mean:             ?[]const u8 = null,
            coefficient:                   ?[]const u8 = null,
            coefficient_mean:              ?[]const u8 = null,
            kinetic_energy:                ?[]const u8 = null,
            kinetic_energy_mean:           ?[]const u8 = null,
            momentum:                      ?[]const u8 = null,
            momentum_mean:                 ?[]const u8 = null,
            population:                    ?[]const u8 = null,
            population_mean:               ?[]const u8 = null,
            position:                      ?[]const u8 = null,
            position_mean:                 ?[]const u8 = null,
            potential_energy:              ?[]const u8 = null,
            potential_energy_mean:         ?[]const u8 = null,
            time_derivative_coupling:      ?[]const u8 = null,
            time_derivative_coupling_mean: ?[]const u8 = null,
            total_energy:                  ?[]const u8 = null,
            total_energy_mean:             ?[]const u8 = null
        };

        derivative_step: T = 0.001,
        iterations: u32,
        seed: u32 = 1,
        time_derivative_coupling: ?[]const u8 = "npi",
        time_step: T,
        trajectories: u32,

        hamiltonian: Hamiltonian, fewest_switches: ?FewestSwitches = null, landau_zener: ?LandauZener = null, spin_mapping: ?SpinMapping = null,

        initial_conditions: InitialConditions, log_intervals: LogIntervals = .{}, write: Write = .{},
    };
}

/// The CI options
pub fn ConfigurationInteractionOptions(comptime T: type) type {
    return struct {
        pub const Gradient = struct {
            numeric: bool = true,
            step: T = 1e-3,
        };
        pub const Hessian = struct {
            numeric: bool = true,
            freq: bool = true,
            step: T = 1e-3,
        };
        pub const Optimize = struct {
            maxiter: u32 = 100,
            threshold: T = 1e-6,
            step: T = 1,
        };
        pub const Print = struct {
            hessian: bool = false,
        };
        pub const Write = struct {
            gradient: ?[]const u8 = null,
            hessian:  ?[]const u8 = null
        };

        active_space: ?[]const usize = null,

        hartree_fock: HartreeFockOptions(T) = .{}, gradient: ?Gradient = null, hessian: ?Hessian = null, optimize: ?Optimize = null, print: Print = .{}, write: Write = .{}
    };
}

/// The Fibonacci options.
pub fn FibonacciOptions(comptime T: type) type {
    return struct {
        count: u32,
        log_interval: u32 = 1,
        output: ?[]const u8 = null,
        start: []const T = &[_]T{0, 1}
    };
}

/// The Hartree-Fock options
pub fn HartreeFockOptions(comptime T: type) type {
    return struct {
        pub const Integral = struct {
            overlap: ?[]const u8 = null,
            kinetic: ?[]const u8 = null,
            nuclear: ?[]const u8 = null,
            coulomb: ?[]const u8 = null,
            basis:   ?[]const u8 = null,
        };
        pub const Gradient = struct {
            numeric: bool = true,
            step: T = 1e-3,
        };
        pub const Hessian = struct {
            numeric: bool = true,
            freq: bool = true,
            step: T = 1e-3,
        };
        pub const Optimize = struct {
            maxiter: u32 = 100,
            threshold: T = 1e-6,
            step: T = 1,
        };
        pub const Print = struct {
            hessian: bool = false,
        };
        pub const System = struct {
            atoms: ?[]const u8 = null,
            charge: i32 = 0,
            coords: ?[]const [3]T = null
        };
        pub const Write = struct {
            coefficient: ?[]const u8 = null,
            density:     ?[]const u8 = null,
            fock:        ?[]const u8 = null,
            nuclear:     ?[]const u8 = null,
            kinetic:     ?[]const u8 = null,
            coulomb:     ?[]const u8 = null,
            overlap:     ?[]const u8 = null,
            gradient:    ?[]const u8 = null,
            hessian:     ?[]const u8 = null
        };

        system_file: ?[]const u8 = null,
        threshold: T = 1e-12,
        maxiter: u32 = 100,
        dsize: ?u32 = 5,
        generalized: bool = false,
        mulliken: bool = true,
        direct: bool = false,

        integral: Integral = .{}, gradient: ?Gradient = null, hessian: ?Hessian = null, optimize: ?Optimize = null, print: Print = .{}, system: System = .{}, write: Write = .{}
    };
}

/// Option struct for the matrix manipulation target.
pub fn MatrixOptions(comptime T: type) type {
    return struct {
        pub const Multiply = struct {
        };
        pub const Random = struct {
            dims: [2]usize,
            distribution: []const u8 = "normal",
            parameters: [2]T = .{0, 1},
            seed: usize = 1,
            count: usize = 1,
        };

        inputs: ?[]const []const u8 = null, outputs: ?[]const []const u8 = null, print: bool = true,

        multiply: ?Multiply = null, random: ?Random = null
    };
}

/// Option struct for the model potential run target
pub fn ModelPotentialOptions(comptime T: type) type {
    return struct {
        pub const Hamiltonian = struct {
            dims: ?u32 = null, matrix: ?[]const []const []const u8 = null, name: ?[]const u8 = null, file: ?[]const u8 = null
        };
        pub const ValuePair = struct {
            index: u32 = 0, value: T = 0
        };

        adiabatic: bool,
        limits: []const T,
        output: []const u8 = "POTENTIAL.mat",
        points: u32,

        constant: []const ValuePair = &[_]ValuePair{}, hamiltonian: Hamiltonian,
    };
}

/// The Moller-Plesset options
pub fn MollerPlessetOptions(comptime T: type) type {
    return struct {
        pub const Gradient = struct {
            numeric: bool = true,
            step: T = 1e-3,
        };
        pub const Hessian = struct {
            numeric: bool = true,
            freq: bool = true,
            step: T = 1e-3,
        };
        pub const Optimize = struct {
            maxiter: u32 = 100,
            threshold: T = 1e-6,
            step: T = 1,
        };
        pub const Print = struct {
            hessian: bool = false,
        };
        pub const Write = struct {
            gradient: ?[]const u8 = null,
            hessian:  ?[]const u8 = null
        };

        order: u32 = 2,

        hartree_fock: HartreeFockOptions(T) = .{}, gradient: ?Gradient = null, hessian: ?Hessian = null, optimize: ?Optimize = null, print: Print = .{}, write: Write = .{}
    };
}

/// The prime generation and checking options
pub fn PrimeOptions(comptime T: type) type {
    return struct {
        pub const Generate = struct {
            count: u32, log_interval: u32 = 1, output: ?[]const u8 = null
        };

        mode: []const u8,
        generate: ?Generate,
        number: T,
    };
}

/// The quantum dynamics options struct
pub fn QuantumDynamicsOptions(comptime T: type) type {
    return struct {
        pub const BohmianDynamics = struct {
            trajectories: usize = 100,
            seed: usize = 1,
        };
        pub const Grid = struct {
            limits: []const T, points: u32
        };
        pub const Hamiltonian = struct {
            dims: ?u32 = null, matrix: ?[]const []const []const u8 = null, cap: ?[]const u8 = null, name: ?[]const u8 = null, file: ?[]const u8 = null
        };
        pub const InitialConditions = struct {
            position: []const T, momentum: []const T, gamma: T = 2, state: u32, mass: T, adiabatic: bool = false
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Spectrum = struct {
            gaussian_window_exponent: T = 0.001,
            nearest_power_of_two: u32 = 2,
            flip: bool = true
        };
        pub const Write = struct {
            autocorrelation_function:             ?[]const u8 = null,
            bloch_vector:                         ?[]const u8 = null,
            bohm_momentum:                        ?[]const u8 = null,
            bohm_momentum_mean:                   ?[]const u8 = null,
            bohm_position:                        ?[]const u8 = null,
            bohm_position_mean:                   ?[]const u8 = null,
            density:                              ?[]const u8 = null,
            kinetic_energy:                       ?[]const u8 = null,
            momentum:                             ?[]const u8 = null,
            population:                           ?[]const u8 = null,
            position:                             ?[]const u8 = null,
            potential_energy:                     ?[]const u8 = null,
            spectrum:                             ?[]const u8 = null,
            total_energy:                         ?[]const u8 = null,
            transformed_autocorrelation_function: ?[]const u8 = null,
            wavefunction:                         ?[]const u8 = null,
        };

        adiabatic: bool,
        iterations: u32,
        mode: []const u32,
        time_step: T,

        bohmian_dynamics: ?BohmianDynamics = null, grid: Grid, hamiltonian: Hamiltonian, initial_conditions: InitialConditions, log_intervals: LogIntervals = .{}, spectrum: Spectrum = .{}, write: Write = .{}
    };
}

/// The sort options
pub fn SortOptions(comptime _: type) type {
    return struct {
        algorithm: []const u8 = "bubble", column: usize = 0, input: []const u8, output: ?[]const u8 = null, print: bool = true
    };
}

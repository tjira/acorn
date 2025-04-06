//! Module where all the output structs are stored.

const std = @import("std");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The classical dynamics options
pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            position_mean: []const T,
            position_std:  []const T,
            momentum_mean: []const T,
            momentum_std:  []const T,
            state:         []const T,
            mass:          []const T
        };
        pub const FewestSwitches = struct {
            quantum_substep: u32 = 10,
            decoherence_alpha: ?T = null
        };
        pub const LandauZener = struct {
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

        adiabatic: bool,
        derivative_step: T = 0.001,
        iterations: u32,
        potential: []const u8,
        seed: u32 = 1,
        time_derivative_coupling: ?[]const u8 = "numeric",
        time_step: T,
        trajectories: u32,

        fewest_switches: ?FewestSwitches = null, landau_zener: ?LandauZener = null, spin_mapping: ?SpinMapping = null,

        initial_conditions: InitialConditions, log_intervals: LogIntervals = .{}, write: Write = .{},
    };
}

/// The CI options
pub fn ConfigurationInteractionOptions(comptime T: type) type {
    return struct {
        excitation: ?[]const u32 = null,

        hartree_fock: HartreeFockOptions(T) = .{},
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
            basis:   ?[]const u8 = "sto-3g"
        };
        pub const System = struct {
            atoms: ?[]const u8 = null,
            charge: u32 = 0,
            coords: ?[]const [3]T = null,
            multiplicity: u32 = 1
        };
        pub const Write = struct {
            coefficient_ao: ?[]const u8 = null,
            density_ao:     ?[]const u8 = null,
            nuclear_ao:     ?[]const u8 = null,
            kinetic_ao:     ?[]const u8 = null,
            coulomb_ao:     ?[]const u8 = null,
            overlap_ao:     ?[]const u8 = null,
        };

        system_file: ?[]const u8 = null,
        threshold: T = 1e-12,
        maxiter: u32 = 100,
        dsize: ?u32 = 5,

        integral: Integral = .{}, system: System = .{}, write: Write = .{}
    };
}

/// Option struct for the model potential run target
pub fn ModelPotentialOptions(comptime T: type) type {
    return struct {
        pub const ValuePair = struct {
            index: u32 = 0, value: T = 0
        };

        adiabatic: bool,
        limits: []const T,
        output: []const u8 = "POTENTIAL.mat",
        points: u32,
        potential: []const u8,

        constant: []const ValuePair = &[_]ValuePair{}
    };
}

/// The Moller-Plesset options
pub fn MollerPlessetOptions(comptime T: type) type {
    return struct {
        order: u32 = 2,

        hartree_fock: HartreeFockOptions(T) = .{},
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
        pub const Grid = struct {
            limits: []const T, points: u32
        };
        pub const InitialConditions = struct {
            position: []const T, momentum: []const T, gamma: T, state: u32, mass: T
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
            transformed_autocorrelation_function: ?[]const u8 = null,
            kinetic_energy:                       ?[]const u8 = null,
            momentum:                             ?[]const u8 = null,
            population:                           ?[]const u8 = null,
            position:                             ?[]const u8 = null,
            potential_energy:                     ?[]const u8 = null,
            spectrum:                             ?[]const u8 = null,
            total_energy:                         ?[]const u8 = null,
            wavefunction:                         ?[]const u8 = null
        };

        adiabatic: bool,
        iterations: u32,
        mode: []const u32,
        potential: []const u8,
        time_step: T,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals = .{}, spectrum: Spectrum = .{}, write: Write = .{}
    };
}

/// The sort options
pub fn SortOptions() type {
    return struct {
        input:     []const u8,
        algorithm: []const u8,
        output:    []const u8
    };
}

//! https://pubs.acs.org/doi/10.1021/jp512221x

const std = @import("std"); const builtin = @import("builtin");

const cdn = @import("acorn").cdn;
const mpt = @import("acorn").mpt;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt_fssh = cdn.ClassicalDynamicsOptions(f64){
        .log_intervals = .{
            .trajectory = 100, .iteration = 500
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = "POPULATION_FSSH.mat",
            .position_mean = null,
            .potential_energy_mean = null,
            .time_derivative_coupling_mean = null,
            .total_energy_mean = null
        },

        .initial_conditions = .{
            .position_mean = &[_]f64{0, 0, 0, 0, 0, 0, 0, 0},
            .position_std  = &[_]f64{11.9739, 11.7264, 9.57297, 8.90160, 8.82918, 8.09775, 7.89327, 7.82045},
            .momentum_mean = &[_]f64{0, 0, 0, 0, 0, 0, 0, 0},
            .momentum_std  = &[_]f64{0.0400449, 0.0411942, 0.0520359, 0.0560927, 0.0565616, 0.0617271, 0.0633333, 0.0639249},
            .state = 2, .mass = 1
        },
        .fewest_switches = .{
            .quantum_substep = 10, .decoherence_alpha = null
        },

        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3500,
        .potential = "uracil1D_1",
        .seed = 1,
        .time_derivative_coupling = "numeric",
        .time_step = 0.012,
        .trajectories = 1000,
    };

    const opt_potential_10 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q10.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 1, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 6, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_12 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q12.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 6, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_18 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q18.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 1, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 6, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_20 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q20.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 1, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 6, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_21 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q21.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 1, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 6, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_24 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q24.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 1, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 6, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_25 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q25.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 1, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 7, .value = 0}
        }
    };

    const opt_potential_26 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4, 4},
        .output = "POTENTIAL_Q26.mat",
        .points = 1024,
        .potential = "uracil1D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0, .value = 0},
            .{.index = 1, .value = 0},
            .{.index = 2, .value = 0},
            .{.index = 3, .value = 0},
            .{.index = 4, .value = 0},
            .{.index = 5, .value = 0},
            .{.index = 6, .value = 0}
        }
    };

    try mpt.write(f64, opt_potential_10, allocator);
    try mpt.write(f64, opt_potential_12, allocator);
    try mpt.write(f64, opt_potential_18, allocator);
    try mpt.write(f64, opt_potential_20, allocator);
    try mpt.write(f64, opt_potential_21, allocator);
    try mpt.write(f64, opt_potential_24, allocator);
    try mpt.write(f64, opt_potential_25, allocator);
    try mpt.write(f64, opt_potential_26, allocator);

    const output_fssh = try cdn.run(f64, opt_fssh , true, allocator); defer output_fssh.deinit();
}

//! https://pubs.acs.org/doi/10.1021/jp512221x
//! https://arxiv.org/abs/2104.06638

const std = @import("std"); const builtin = @import("builtin");

const cdn = @import("acorn").cdn;
const mpt = @import("acorn").mpt;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const omega = &[_]f64{388, 560, 734, 770, 771, 1193, 1228, 1383, 1406, 1673, 1761, 1794};

    var std_x: [12]f64 = undefined; var std_p: [12]f64 = undefined; for (0..omega.len) |i| {
        const a = std.math.sqrt(219474.6305 * 1.0 / 1.0 / omega[i]); std_x[i] = a / std.math.sqrt2; std_p[i] = 1 / a / std.math.sqrt2;
    }

    const opt_fssh_12D = cdn.ClassicalDynamicsOptions(f64){
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
            .position_mean = &[_]f64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            .position_std  = &std_x,
            .momentum_mean = &[_]f64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            .momentum_std  = &std_p,
            .state = 2, .mass = 1
        },
        .fewest_switches = .{
            .quantum_substep = 100, .decoherence_alpha = null
        },

        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3500,
        .potential = "uracil12D_1",
        .seed = 1,
        .time_derivative_coupling = "numeric",
        .time_step = 0.012,
        .trajectories = 1000,
    };

    const opt_potential_adia_3 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_Q3.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_7 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3, 3},
        .output = "ADIABATIC_POTENTIAL_Q7.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_10 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-2, 2},
        .output = "ADIABATIC_POTENTIAL_Q10.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value =  0},
            .{.index = 1,  .value =  0},
            .{.index = 3,  .value =  0},
            .{.index = 4,  .value =  0},
            .{.index = 5,  .value =  0},
            .{.index = 6,  .value =  0},
            .{.index = 7,  .value =  0},
            .{.index = 8,  .value = -4},
            .{.index = 9,  .value =  0},
            .{.index = 10, .value =  0},
            .{.index = 11, .value =  0}
        }
    };

    const opt_potential_adia_11 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3, 3},
        .output = "ADIABATIC_POTENTIAL_Q11.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_12 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-2, 2},
        .output = "ADIABATIC_POTENTIAL_Q12.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value =  0},
            .{.index = 1,  .value =  0},
            .{.index = 2,  .value =  0},
            .{.index = 3,  .value =  0},
            .{.index = 5,  .value =  0},
            .{.index = 6,  .value =  0},
            .{.index = 7,  .value =  0},
            .{.index = 8,  .value = -4},
            .{.index = 9,  .value =  0},
            .{.index = 10, .value =  0},
            .{.index = 11, .value =  0}
        }
    };

    const opt_potential_adia_18 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-2, 2},
        .output = "ADIABATIC_POTENTIAL_Q18.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_19 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3.5, 3.5},
        .output = "ADIABATIC_POTENTIAL_Q19.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_20 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3.5, 3.5},
        .output = "ADIABATIC_POTENTIAL_Q20.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_21 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_Q21.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_24 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3.5, 3.5},
        .output = "ADIABATIC_POTENTIAL_Q24.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 10, .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_25 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_Q25.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 11, .value = 0}
        }
    };

    const opt_potential_adia_26 = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_Q26.mat",
        .points = 1024,
        .potential = "uracil12D_1",

        .constant = &[_]mpt.ModelPotentialOptions(f64).ValuePair{
            .{.index = 0,  .value = 0},
            .{.index = 1,  .value = 0},
            .{.index = 2,  .value = 0},
            .{.index = 3,  .value = 0},
            .{.index = 4,  .value = 0},
            .{.index = 5,  .value = 0},
            .{.index = 6,  .value = 0},
            .{.index = 7,  .value = 0},
            .{.index = 8,  .value = 0},
            .{.index = 9,  .value = 0},
            .{.index = 10, .value = 0}
        }
    };

    try mpt.write(f64, opt_potential_adia_3,  allocator);
    try mpt.write(f64, opt_potential_adia_7,  allocator);
    try mpt.write(f64, opt_potential_adia_10, allocator);
    try mpt.write(f64, opt_potential_adia_11, allocator);
    try mpt.write(f64, opt_potential_adia_12, allocator);
    try mpt.write(f64, opt_potential_adia_18, allocator);
    try mpt.write(f64, opt_potential_adia_19, allocator);
    try mpt.write(f64, opt_potential_adia_20, allocator);
    try mpt.write(f64, opt_potential_adia_21, allocator);
    try mpt.write(f64, opt_potential_adia_24, allocator);
    try mpt.write(f64, opt_potential_adia_25, allocator);
    try mpt.write(f64, opt_potential_adia_26, allocator);

    const output_fssh_12D = try cdn.run(f64, opt_fssh_12D , true, allocator); defer output_fssh_12D.deinit();
}

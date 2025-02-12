//! https://pubs.acs.org/doi/10.1021/jp512221x
//! https://pubs.acs.org/doi/10.1021/acs.jpclett.1c04132
//! https://arxiv.org/abs/2104.06638

const std = @import("std"); const builtin = @import("builtin");

const classical_dynamics = @import("acorn").classical_dynamics;
const constant           = @import("acorn").constant          ;
const model_potential    = @import("acorn").model_potential   ;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    var omega_8D:  [ 8]f64 = .{          734,      771, 1193,       1383, 1406, 1673, 1761, 1794}; for (0..omega_8D.len ) |i| omega_8D[i]  /= constant.EV2RCM * constant.AU2EV;
    var omega_12D: [12]f64 = .{388, 560, 734, 770, 771, 1193, 1228, 1383, 1406, 1673, 1761, 1794}; for (0..omega_12D.len) |i| omega_12D[i] /= constant.EV2RCM * constant.AU2EV;

    var std_x_8D: [8]f64 = undefined; var std_p_8D: [8]f64 = undefined; for (0..omega_8D.len) |i| {
        std_x_8D[i] = std.math.sqrt(0.5 / omega_8D[i]); std_p_8D[i] = 0.08;
    }

    const opt_fssh_8D = classical_dynamics.ClassicalDynamicsOptions(f64){
        .log_intervals = .{
            .trajectory = 100, .iteration = 50
        },
        .write = .{
            .population_mean = "POPULATION_FSSH_8D.mat"
        },

        .initial_conditions = .{
            .position_mean = &[_]f64{0, 0, 0, 0, 0, 0, 0, 0},
            .position_std  = &std_x_8D,
            .momentum_mean = &[_]f64{0, 0, 0, 0, 0, 0, 0, 0},
            .momentum_std  = &std_p_8D,
            .state = &[_]f64{0.01, 0.05, 0.94, 0}, .mass = &[_]f64{1, 1, 1, 1, 1, 1, 1, 1}
        },
        .fewest_switches = .{
            .quantum_substep = 100, .decoherence_alpha = null
        },

        .adiabatic = true,
        .iterations = 250,
        .potential = "uracil8D_1",
        .time_step = 10,
        .trajectories = 1000
    };

    var opt_lzsh_8D = opt_fssh_8D; opt_lzsh_8D.fewest_switches = null; opt_lzsh_8D.landau_zener = .{}; opt_lzsh_8D.write.population_mean = "POPULATION_LZSH_8D.mat";
    var opt_mash_8D = opt_fssh_8D; opt_mash_8D.fewest_switches = null; opt_mash_8D.spin_mapping = .{}; opt_mash_8D.write.population_mean = "POPULATION_MASH_8D.mat";

    var opt_ktsh_8D = opt_fssh_8D; opt_ktsh_8D.time_derivative_coupling = "baeckan"; opt_ktsh_8D.write.population_mean = "POPULATION_KTSH_8D.mat";

    const opt_potential_adia_3 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q3.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_7 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3, 3},
        .output = "ADIABATIC_POTENTIAL_12D_Q7.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_10 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-2, 2},
        .output = "ADIABATIC_POTENTIAL_12D_Q10.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_11 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3, 3},
        .output = "ADIABATIC_POTENTIAL_12D_Q11.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_12 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-2, 2},
        .output = "ADIABATIC_POTENTIAL_12D_Q12.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_18 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-2, 2},
        .output = "ADIABATIC_POTENTIAL_12D_Q18.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_19 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3.5, 3.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q19.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_20 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3.5, 3.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q20.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_21 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q21.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_24 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-3.5, 3.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q24.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_25 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q25.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    const opt_potential_adia_26 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-4.5, 4.5},
        .output = "ADIABATIC_POTENTIAL_12D_Q26.mat",
        .points = 1024,
        .potential = "uracilDimless12D_1",

        .constant = &[_]model_potential.ModelPotentialOptions(f64).ValuePair{
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

    try model_potential.write(f64, opt_potential_adia_3,  allocator);
    try model_potential.write(f64, opt_potential_adia_7,  allocator);
    try model_potential.write(f64, opt_potential_adia_10, allocator);
    try model_potential.write(f64, opt_potential_adia_11, allocator);
    try model_potential.write(f64, opt_potential_adia_12, allocator);
    try model_potential.write(f64, opt_potential_adia_18, allocator);
    try model_potential.write(f64, opt_potential_adia_19, allocator);
    try model_potential.write(f64, opt_potential_adia_20, allocator);
    try model_potential.write(f64, opt_potential_adia_21, allocator);
    try model_potential.write(f64, opt_potential_adia_24, allocator);
    try model_potential.write(f64, opt_potential_adia_25, allocator);
    try model_potential.write(f64, opt_potential_adia_26, allocator);

    const output_fssh_8D = try classical_dynamics.run(f64, opt_fssh_8D , true, allocator); defer output_fssh_8D.deinit();
    const output_lzsh_8D = try classical_dynamics.run(f64, opt_lzsh_8D , true, allocator); defer output_lzsh_8D.deinit();
    const output_ktsh_8D = try classical_dynamics.run(f64, opt_ktsh_8D , true, allocator); defer output_ktsh_8D.deinit();
    const output_mash_8D = try classical_dynamics.run(f64, opt_mash_8D , true, allocator); defer output_mash_8D.deinit();
}

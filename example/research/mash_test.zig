const std = @import("std"); const builtin = @import("builtin");

const classical_dynamics = @import("acorn").classical_dynamics;
const quantum_dynamics   = @import("acorn").quantum_dynamics  ;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt_exact_tully1D_2 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 50
        },

        .grid = .{
            .limits = &[_]f64{-48, 144}, .points = 4096
        },

        .initial_conditions = .{
            .position = &[_]f64{-15}, .momentum = &[_]f64{35}, .gamma = 0.5, .state = 0, .mass = 2000
        },

        .write = .{
            .population = "POPULATION_EXACT_tully1D_2.mat",
        },

        .adiabatic = true,
        .iterations = 207,
        .mode = &[_]u32{0, 1},
        .potential = "tully1D_2",
        .time_step = 10,
    };

    const opt_exact_tully1D_3 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 50
        },

        .grid = .{
            .limits = &[_]f64{-48, 144}, .points = 4096
        },

        .initial_conditions = .{
            .position = &[_]f64{-15}, .momentum = &[_]f64{10}, .gamma = 0.5, .state = 1, .mass = 2000
        },

        .write = .{
            .population = "POPULATION_EXACT_tully1D_3.mat",
        },

        .adiabatic = true,
        .iterations = 621,
        .mode = &[_]u32{0, 1},
        .potential = "tully1D_3",
        .time_step = 10,
    };

    const opt_mash_tully1D_2 = classical_dynamics.ClassicalDynamicsOptions(f64){
        .log_intervals = .{
            .trajectory = 100, .iteration = 500
        },

        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{1.0},
            .momentum_mean = &[_]f64{35.0},
            .momentum_std  = &[_]f64{0.5},
            .state = &[_]f64{1, 0}, .mass = &[_]f64{2000}
        },

        .spin_mapping = .{
            .resample = null
        },

        .write = .{
            .population_mean = "POPULATION_MASH_tully1D_2.mat",
        },

        .adiabatic = true,
        .iterations = 2070,
        .potential = "tully1D_2",
        .time_step = 1,
        .trajectories = 10000,
    };

    const opt_mash_tully1D_3 = classical_dynamics.ClassicalDynamicsOptions(f64){
        .log_intervals = .{
            .trajectory = 100, .iteration = 500
        },

        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{1.0},
            .momentum_mean = &[_]f64{10.0},
            .momentum_std  = &[_]f64{0.5},
            .state = &[_]f64{1, 0}, .mass = &[_]f64{2000}
        },

        .spin_mapping = .{
            .resample = .{.reflect = true}
        },

        .write = .{
            .population_mean = "POPULATION_MASH_tully1D_3.mat",
        },

        .adiabatic = true,
        .iterations = 6210,
        .potential = "tully1D_3",
        .time_step = 1,
        .trajectories = 10000,
    };

    const output_exact_tully1D_2 =   try quantum_dynamics.run(f64, opt_exact_tully1D_2, true, allocator); defer output_exact_tully1D_2.deinit();
    const output_exact_tully1D_3 =   try quantum_dynamics.run(f64, opt_exact_tully1D_3, true, allocator); defer output_exact_tully1D_3.deinit();
    const output_mash_tully1D_2  = try classical_dynamics.run(f64, opt_mash_tully1D_2 , true, allocator); defer  output_mash_tully1D_2.deinit();
    const output_mash_tully1D_3  = try classical_dynamics.run(f64, opt_mash_tully1D_3 , true, allocator); defer  output_mash_tully1D_3.deinit();
}

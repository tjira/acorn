const std = @import("std"); const builtin = @import("builtin");

const classical_dynamics = @import("acorn").classical_dynamics;
const quantum_dynamics   = @import("acorn").quantum_dynamics  ;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt_exact = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .population = "POPULATION_EXACT.mat",
        },

        .grid = .{
            .limits = &[_]f64{-48, 144}, .points = 4096
        },
        .initial_conditions = .{
            .position = &[_]f64{-15}, .momentum = &[_]f64{25}, .gamma = 2, .state = 0, .mass = 2000
        },

        .adiabatic = true,
        .iterations = 350,
        .mode = &[_]u32{0, 1},
        .potential = "tully1D_3",
        .time_step = 10,
    };

    const opt_fssh = classical_dynamics.ClassicalDynamicsOptions(f64){
        .log_intervals = .{
            .trajectory = 100, .iteration = 500
        },
        .write = .{
            .population_mean = "POPULATION_FSSH.mat",
        },

        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{0.5},
            .momentum_mean = &[_]f64{25.0},
            .momentum_std  = &[_]f64{1.0},
            .state = 1, .mass = 2000
        },
        .fewest_switches = .{
            .quantum_substep = 10, .decoherence_alpha = null
        },

        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3500,
        .potential = "tully1D_3",
        .seed = 1,
        .time_derivative_coupling = "numeric",
        .time_step = 1,
        .trajectories = 1000,
    };

    var opt_mash = opt_fssh;
    opt_mash.write.population_mean = "POPULATION_MASH.mat";
    opt_mash.fewest_switches = null;
    opt_mash.spin_mapping = .{
        .quantum_substep = 10, .sample_bloch_vector = false,
    };

    const output_exact =   try quantum_dynamics.run(f64, opt_exact, true, allocator); defer output_exact.deinit();
    const output_fssh  = try classical_dynamics.run(f64, opt_fssh , true, allocator); defer  output_fssh.deinit();
    const output_mash  = try classical_dynamics.run(f64, opt_mash , true, allocator); defer  output_mash.deinit();
}

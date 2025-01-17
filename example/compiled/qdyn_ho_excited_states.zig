const std = @import("std"); const builtin = @import("builtin");

const qdn = @import("acorn").qdn;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt = qdn.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },
        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },
        .log_intervals = .{
            .iteration = 2000
        },
        .write = .{
            .autocorrelation_function = null,
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null,
            .wavefunction = null
        },

        .adiabatic = false,
        .iterations = 10000,
        .mode = &[_]u32{8, 0},
        .potential = "harmonic1D_1",
        .time_step = 0.01,
    };

    const output = try qdn.run(f64, opt, true, allocator); defer output.deinit();
}

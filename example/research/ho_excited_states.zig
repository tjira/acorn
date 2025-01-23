const std = @import("std"); const builtin = @import("builtin");

const qdn = @import("acorn").qdn;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt_imaginary = qdn.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 2000
        },
        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },
        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },

        .adiabatic = false,
        .iterations = 10000,
        .mode = &[_]u32{8, 0},
        .potential = "harmonic1D_1",
        .time_step = 0.01,
    };

    var opt_real = opt_imaginary; opt_real.mode = &[_]u32{0, 1};

    opt_real.write.spectrum                           = "SPECTRUM.mat";
    opt_real.write.autocorrelation_function           =      "ACF.mat";

    const output_imaginary = try qdn.run(f64, opt_imaginary, true, allocator); defer output_imaginary.deinit();
    const output_real      = try qdn.run(f64, opt_real,      true, allocator); defer      output_real.deinit();
}

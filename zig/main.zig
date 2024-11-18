const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );
const qdn = @import("quantumdynamics.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    var timer = try std.time.Timer.start();

    // const opt = qdn.QuantumDynamicsOptions(f64){
    //     .iterations = 350,
    //     .time_step = 10,
    //     .imaginary = false,
    //     .grid = .{
    //         .limits = &[_]f64{-16, 24},
    //         .points = 256
    //     },
    //     .initial_conditions = .{
    //         .position = &[_]f64{-10},
    //         .momentum = &[_]f64{15},
    //         .state = 1,
    //         .mass = 2000
    //     },
    //     .log_intervals = .{
    //         .iteration = 50
    //     },
    //     .write = .{
    //         .population = null
    //     },
    //     .potential = mpt.doubleState1D_1,
    // };
    const opt = qdn.QuantumDynamicsOptions(f64){
        .iterations = 100,
        .time_step = 0.1,
        .imaginary = true,
        .grid = .{
            .limits = &[_]f64{-8, 8},
            .points = 64
        },
        .initial_conditions = .{
            .position = &[_]f64{1},
            .momentum = &[_]f64{0},
            .state = 0,
            .mass = 1
        },
        .log_intervals = .{
            .iteration = 10
        },
        .write = .{
            .population = null
        },
        .potential = mpt.harmonic1D_1,
    };
    try qdn.run(f64, opt, allocator);

    const pot_start: f64 = -16; const pot_end: f64 = 16; const pot_points: f64 = 1024;
    try mpt.write(f64, "POTENTIAL-DIABATIC.mat",  opt.potential, pot_start, pot_end, pot_points, false, allocator);
    try mpt.write(f64, "POTENTIAL-ADIABATIC.mat", opt.potential, pot_start, pot_end, pot_points, true,  allocator);

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});
}

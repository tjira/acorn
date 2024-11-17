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

    const pot_start: f64 = -16; const pot_end: f64 = 16; const pot_points: f64 = 1024;

    const qdyn_opt = qdn.QuantumDynamicsOptions(f64){
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
            .iteration = 500
        },
        .potential = mpt.harmonic1D_1,
    };

    // const qdyn_opt = qdn.QuantumDynamicsOptions(f64){
    //     .iterations = 350,
    //     .time_step = 10,
    //     .imaginary = false,
    //     .grid = .{
    //         .limits = &[_]f64{-16, 16},
    //         .points = 256
    //     },
    //     .initial_conditions = .{
    //         .position = &[_]f64{-10},
    //         .momentum = &[_]f64{15},
    //         .state = 1,
    //         .mass = 2000
    //     },
    //     .log_intervals = .{
    //         .iteration = 500
    //     },
    //     .potential = mpt.doubleState1D_1,
    // };

    const cdyn_opt = cdn.ClassicalDynamicsOptions(f64){
        .adiabatic = false,
        .iterations = 3500,
        .seed = 1,
        .time_step = 1,
        .derivative_step = 0.001,
        .trajectories = 1000,
        .initial_conditions = .{
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .state = 1,
            .mass = 2000
        },
        .log_intervals = .{
            .trajectory = 100,
            .iteration = 500
        },
        .potential = mpt.tripleState1D_1,
    };

    try mpt.write(f64, "POTENTIAL-DIABATIC.mat",  cdyn_opt.potential, pot_start, pot_end, pot_points, false, allocator);
    try mpt.write(f64, "POTENTIAL-ADIABATIC.mat", cdyn_opt.potential, pot_start, pot_end, pot_points, true,  allocator);

    // try cdn.run(f64, cdyn_opt, allocator);
    try qdn.run(f64, qdyn_opt, allocator);

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});

    // const a = try Matrix(std.math.Complex(f64)).init(6, 1, allocator);
    // const b = try Matrix(std.math.Complex(f64)).init(6, 1, allocator);
    //
    // a.ptr(0, 0).* = std.math.Complex(f64).init(1, 0);
    // a.ptr(1, 0).* = std.math.Complex(f64).init(8, 0);
    // a.ptr(2, 0).* = std.math.Complex(f64).init(3, 0);
    // a.ptr(3, 0).* = std.math.Complex(f64).init(0, 0);
    // a.ptr(4, 0).* = std.math.Complex(f64).init(1, 0);
    // a.ptr(5, 0).* = std.math.Complex(f64).init(2, 0);
    //
    // ftr.dft(f64, b.data, a.data, &[_]usize{a.rows}, -1);
    //
    // for (a.data) |e| std.debug.print("{}\n", .{e});
    // std.debug.print("\n", .{});
    // for (b.data) |e| std.debug.print("{}\n", .{e});
    // std.debug.print("\n", .{});
    //
    // ftr.dft(f64, a.data, b.data, &[_]usize{a.rows}, 1);
    // for (a.data) |e| std.debug.print("{}\n", .{e});
}

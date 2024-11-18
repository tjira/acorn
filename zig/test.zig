const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const mat = @import("matrix.zig"           );
const mpt = @import("modelpotential.zig"   );
const qdn = @import("quantumdynamics.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.testing.allocator; const tol = 1e-14;

test "cdnDoubleState1D_1-diabatic" {
    const opt = cdn.ClassicalDynamicsOptions(f64){
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
        .write = .{
            .population = null
        },
        .potential = mpt.doubleState1D_1,
    };
    try cdn.run(f64, opt, allocator);
    try std.testing.expect(true);
}

test "qdnDoubleState1D_1-diabatic-real" {
    const opt = qdn.QuantumDynamicsOptions(f64){
        .iterations = 350,
        .time_step = 10,
        .imaginary = false,
        .grid = .{
            .limits = &[_]f64{-16, 24},
            .points = 256
        },
        .initial_conditions = .{
            .position = &[_]f64{-10},
            .momentum = &[_]f64{15},
            .state = 1,
            .mass = 2000
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .population = null
        },
        .potential = mpt.doubleState1D_1,
    };
    try qdn.run(f64, opt, allocator);
    try std.testing.expect(true);
}

test "qdnHarmonic1D_1-diabatic-imag" {
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
    try std.testing.expect(true);
}

test "qdnTripleState1D_1-diabatic-real" {
    const opt = qdn.QuantumDynamicsOptions(f64){
        .iterations = 350,
        .time_step = 10,
        .imaginary = false,
        .grid = .{
            .limits = &[_]f64{-16, 24},
            .points = 512
        },
        .initial_conditions = .{
            .position = &[_]f64{-10},
            .momentum = &[_]f64{15},
            .state = 1,
            .mass = 2000
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .population = null
        },
        .potential = mpt.tripleState1D_1,
    };
    try qdn.run(f64, opt, allocator);
    try std.testing.expect(true);
}

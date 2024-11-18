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

    // var map = std.StringHashMap(mpt.PotentialType(f64)).init(allocator);
    //
    // try map.put("tripleState1D_1", mpt.tripleState1D_1);

    // map = map;

    const cdyn_opt = cdn.ClassicalDynamicsOptions(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 1000,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .population = "POPULATION-LZSH.mat"
        },
        .potential = mpt.tripleState1D_1,
    };

    const qdyn_opt = qdn.QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 1024
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .population = "POPULATION-EXACT.mat"
        },
        .potential = mpt.tripleState1D_1,
    };

    try cdn.run(f64, cdyn_opt, allocator); try qdn.run(f64, qdyn_opt, allocator);

    // const opt = qdn.QuantumDynamicsOptions(f64){
    //     .iterations = 100,
    //     .time_step = 0.1,
    //     .imaginary = true,
    //     .grid = .{
    //         .limits = &[_]f64{-8, 8},
    //         .points = 64
    //     },
    //     .initial_conditions = .{
    //         .position = &[_]f64{1},
    //         .momentum = &[_]f64{0},
    //         .state = 0,
    //         .mass = 1
    //     },
    //     .log_intervals = .{
    //         .iteration = 10
    //     },
    //     .write = .{
    //         .population = null
    //     },
    //     .potential = mpt.harmonic1D_1,
    // };

    const pot_start: f64 = -16; const pot_end: f64 = 16; const pot_points: f64 = 1024;
    try mpt.write(f64, "POTENTIAL-DIABATIC.mat",  cdyn_opt.potential, pot_start, pot_end, pot_points, false, allocator);
    try mpt.write(f64, "POTENTIAL-ADIABATIC.mat", cdyn_opt.potential, pot_start, pot_end, pot_points, true,  allocator);

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});
}

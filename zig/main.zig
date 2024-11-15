const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const mpt = @import("modelpotential.zig"   );
const mat = @import("matrix.zig"           );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    var timer = try std.time.Timer.start();

    const cdyn_opt = cdn.ClassicalDynamicsOptions(f64){
        .adiabatic = false,
        .iterations = 3500,
        .seed = 1,
        .time_step = 1,
        .derivative_step = 0.001,
        .trajectories = 1,
        .ic = .{
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .state = 1,
            .mass = 2000
        },
        .li = .{
            .trajectory = 10,
            .iteration = 500
        },
        .potential = mpt.doubleState1D_1,
    };

    try cdn.run(f64, cdyn_opt, allocator);

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});

    var A = try Matrix(f64).init(5, 2, allocator); defer A.deinit(); A.linspace(-4, 4);
    mat.print(f64, A);
}

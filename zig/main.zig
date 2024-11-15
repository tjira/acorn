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
        .trajectories = 1000,
        .ic = .{
            .position_mean = &[_]f64{-15},
            .position_std = &[_]f64{0.5},
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .state = 1,
            .mass = 2000
        },
        .li = .{
            .trajectory = 100,
            .iteration = 500
        },
        .potential = mpt.tripleState1D_1,
    };

    try mpt.write(f64, "POTENTIAL.mat", cdyn_opt.potential, -16, 16, 1024, cdyn_opt.adiabatic, allocator); try cdn.run(f64, cdyn_opt, allocator);

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});


    // const r = try mpt.grid(f64, -16, 16, 1024, 1, allocator); defer r.deinit();
    // const U = try mpt.evaluate(f64, mpt.tripleState1D_1, r, allocator); defer U.deinit();
    //
    // var Ur = try Matrix(f64).init(U.rows, r.cols + U.cols, allocator);
    //
    // r.hjoin(&Ur, U); try Ur.write("POTENTIAL.mat");
    //
    // try Ur.print(std.io.getStdOut().writer());

    // const start = -16; const end = 16; const points = 1024;
    // var Ur = try Matrix(f64).init(points, 4, allocator); defer Ur.deinit();
    // var U = try Matrix(f64).init(3, 3, allocator); defer U.deinit();
    // var r = try Vector(f64).init(1, allocator); defer r.deinit();
    // for (0..points) |i| {
    //     r.ptr(0).* = start + (end - start) / @as(f64, @floatFromInt(points - 1)) * @as(f64, @floatFromInt(i));
    //     Ur.ptr(i, 0).* = r.at(0); mpt.tripleState1D_1(f64, &U, r);
    //     for (0..U.rows) |j| Ur.ptr(i, j + 1).* = U.at(j, j);
    // }
    // try Ur.write("POTENTIAL.mat");

    // var A = try Matrix(f64).init(2, 2, allocator); defer A.deinit();
    // var B = try Matrix(f64).init(2, 2, allocator); defer B.deinit();
    // var C = try Matrix(f64).init(2, 2, allocator); defer C.deinit();
    // A.set(&[_]f64{1, 2, 2, 1}); try mat.eigh(f64, &B, &C, A, 1e-12);
    //
    // try A.print(std.io.getStdOut().writer());
    // try B.print(std.io.getStdOut().writer());
    // try C.print(std.io.getStdOut().writer());
}

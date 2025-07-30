//! Functions to optimize a system with quantum mechanical gradients.

const std = @import("std");

const edf = @import("energy_derivative.zig");
const mat = @import("matrix.zig"           );

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;

/// General function to optimize a system with quantum mechanical gradients.
pub fn optimize(comptime T: type, opt: anytype, system: System(T), efunc: anytype, method: []const u8, print: bool, allocator: std.mem.Allocator) !System(T) {
    if (opt.gradient == null) return error.GradientNotDefined;

    var optsystem = try system.clone();

    if (print) try std.io.getStdOut().writer().print("\n{s} GEOMETRY OPTIMIZATION:\n{s:4} {s:20} {s:4}\n", .{method, "ITER", "GRADIENT NORM", "TIME"});

    for (0..opt.optimize.?.maxiter) |i| {

        if (i == opt.optimize.?.maxiter) return error.MaxIterationsExceeded;

        var timer = try std.time.Timer.start();

        if (print) try std.io.getStdOut().writer().print("{d:4} ", .{i + 1});

        var G = try edf.gradient(T, opt, optsystem, efunc, method, false, allocator);

        if (print) try std.io.getStdOut().writer().print("{d:20.14} {s}\n", .{G.vector().norm(), std.fmt.fmtDuration(timer.read())});

        if (G.vector().norm() < opt.optimize.?.threshold) break;

        mat.muls(T, &G, G, opt.optimize.?.step);

        mat.sub(T, &optsystem.coords, optsystem.coords, G);
    }

    return optsystem;
}

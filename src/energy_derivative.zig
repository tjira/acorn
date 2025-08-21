//! File with functions related to the calculation of energy derivatives from quantum mechanical calculations.

const std = @import("std");

const hlp = @import("helper.zig");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;

/// Function that returns a gradient given the system and energy function.
pub fn gradient(comptime T: type, opt: anytype, system: System(T), efunc: anytype, method: []const u8, print: bool, allocator: std.mem.Allocator) !Matrix(T) {
    var G = try Matrix(T).init(system.coords.rows, 3, allocator);

    if (print) try hlp.print(std.fs.File.stdout(), "\n{s} NUMERICAL GRADIENT:\n{s:7} {s:20} {s:4}\n", .{method, "INDEX", "GRADIENT ELEMENT", "TIME"});

    for (0..system.coords.rows) |i| for (0..3) |j| {

        var timer = try std.time.Timer.start();

        if (print) try hlp.print(std.fs.File.stdout(), "{d:3}/{d:3} ", .{3 * i + j + 1, 3 * system.coords.rows});

        const sysp1 = try system.clone(); defer sysp1.deinit(); sysp1.coords.ptr(i, j).* += opt.gradient.?.step;
        const sysm1 = try system.clone(); defer sysm1.deinit(); sysm1.coords.ptr(i, j).* -= opt.gradient.?.step;

        var outp1 = try efunc(T, opt, sysp1, false, allocator); const Ep1 = outp1.E; outp1.deinit();
        var outm1 = try efunc(T, opt, sysm1, false, allocator); const Em1 = outm1.E; outm1.deinit();

        G.ptr(i, j).* = (Ep1 - Em1) / (2 * opt.gradient.?.step);

        if (print) try hlp.print(std.fs.File.stdout(), "{d:20.14} {D}\n", .{G.at(i, j), timer.read()});
    };

    return G;
}

/// Returns a hessian given the system and energy function.
pub fn hessian(comptime T: type, opt: anytype, system: System(T), efunc: anytype, method: []const u8, print: bool, allocator: std.mem.Allocator) !Matrix(T) {
    var H = try Matrix(T).init(system.coords.rows * system.coords.rows, system.coords.rows * system.coords.rows, allocator); var k: usize = 0;

    if (print) try hlp.print(std.fs.File.stdout(), "\n{s} NUMERICAL HESSIAN:\n{s:7} {s:20} {s:4}\n", .{method, "INDEX", "HESSIAN ELEMENT", "TIME"});

    for (0..H.rows) |i| for (i..H.cols) |j| {

        var timer = try std.time.Timer.start();

        if (print) try hlp.print(std.fs.File.stdout(), "{d:3}/{d:3} ", .{k + 1, H.rows * (H.rows - 1) / 2 + H.rows});

        const sysp2   = try system.clone(); defer   sysp2.deinit(); sysp2.coords.ptr(  i / 3, i % 3).* += opt.hessian.?.step; sysp2.coords.ptr(  j / 3, j % 3).* += opt.hessian.?.step;
        const sysp1m1 = try system.clone(); defer sysp1m1.deinit(); sysp1m1.coords.ptr(i / 3, i % 3).* += opt.hessian.?.step; sysp1m1.coords.ptr(j / 3, j % 3).* -= opt.hessian.?.step;
        const sysm1p1 = try system.clone(); defer sysm1p1.deinit(); sysm1p1.coords.ptr(i / 3, i % 3).* -= opt.hessian.?.step; sysm1p1.coords.ptr(j / 3, j % 3).* += opt.hessian.?.step;
        const sysm2   = try system.clone(); defer   sysm2.deinit(); sysm2.coords.ptr(  i / 3, i % 3).* -= opt.hessian.?.step; sysm2.coords.ptr(  j / 3, j % 3).* -= opt.hessian.?.step;

        var outp2   = try efunc(T, opt, sysp2,   false, allocator); const Ep2   = outp2.E;     outp2.deinit();
        var outp1m1 = try efunc(T, opt, sysp1m1, false, allocator); const Ep1m1 = outp1m1.E; outp1m1.deinit();
        var outm1p1 = try efunc(T, opt, sysm1p1, false, allocator); const Em1p1 = outm1p1.E; outm1p1.deinit();
        var outm2   = try efunc(T, opt, sysm2,   false, allocator); const Em2   = outm2.E;     outm2.deinit();

        H.ptr(i, j).* = (Ep2 - Ep1m1 - Em1p1 + Em2) / (4 * opt.hessian.?.step * opt.hessian.?.step); H.ptr(j, i).* = H.at(i, j); k += 1;

        if (print) try hlp.print(std.fs.File.stdout(), "{d:20.14} {D}\n", .{H.at(i, j), timer.read()});
    };

    return H;
}

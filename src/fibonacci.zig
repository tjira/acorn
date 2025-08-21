//! Module for generating Fibonacci numbers.

const std = @import("std");

const hlp = @import("helper.zig");
const inp = @import("input.zig" );

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The main function for generating Fibonacci numbers.
pub fn run(comptime T: type, opt: inp.FibonacciOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    var timer = try std.time.Timer.start(); var F: Matrix(T) = undefined; var fibs: [2]T = .{opt.start[0], opt.start[1]};

    if (opt.output != null) {
        F = try Matrix(T).init(opt.count, 2, allocator); F.ptr(0, 0).* = 1; F.ptr(0, 1).* = fibs[0]; F.ptr(1, 0).* = 2; F.ptr(1, 1).* = fibs[1];
    }

    if (print) {
        try hlp.print(std.fs.File.stdout(), "\nFIBONACCI NUMBERS\n", .{});
    }

    if (print and opt.count >= 1 and opt.log_interval >= 1) try hlp.print(std.fs.File.stdout(), "0: {d} ({D})\n", .{fibs[0], timer.read()});
    if (print and opt.count >= 2 and opt.log_interval <= 2) try hlp.print(std.fs.File.stdout(), "1: {d} ({D})\n", .{fibs[1], timer.read()});

    for (2..opt.count + 1) |i| {

        fibs[i % 2] = fibs[0] + fibs[1];

        if (opt.output != null) {
            F.ptr(i, 0).* = @as(T, @intCast(i)) + 1; F.ptr(i, 1).* = fibs[i % 2];
        }

        if (print and (i + 1) % opt.log_interval == 0 or i == 1) {
            try hlp.print(std.fs.File.stdout(), "{d}: {d} ({D})\n", .{i, fibs[i % 2], timer.read()});
        }
    }

    if (opt.output != null) {
        try F.write(opt.output.?); F.deinit();
    }
}

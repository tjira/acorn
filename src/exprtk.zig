//! Wrapper for the Exprtk C++ library.

const exprtk = @cImport(@cInclude("exprtk.h"));

const std = @import("std");

const Vector = @import("vector.zig").Vector;

pub fn compile(string: []const u8, nvars: usize, allocator: std.mem.Allocator) !*anyopaque {
    const buffer = try allocator.dupeZ(u8, string); defer allocator.free(buffer); return exprtk.compile(buffer, @as(i32, @intCast(nvars))).?;
}

pub fn evaluate(expr: *anyopaque, vars: Vector(f64)) f64 {
    return exprtk.evaluate(expr, vars.data.ptr);
}

pub fn deinit(expr: *anyopaque) void {
    exprtk.deinit(expr);
}

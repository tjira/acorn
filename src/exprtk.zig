//! Wrapper for the Exprtk C++ library.

const exprtk = @cImport(@cInclude("exprtk.h"));

const std = @import("std");

const Vector = @import("vector.zig").Vector;

/// Compile the algebrainc expression string with the given number of variables. Variables must be in the form of r1, r2, ..., rn.
pub fn compile(string: []const u8, nvars: usize, allocator: std.mem.Allocator) !*anyopaque {
    const buffer = try allocator.dupeZ(u8, string); defer allocator.free(buffer); return exprtk.compile(buffer, nvars).?;
}

/// Evaluate the compiled expression at the specified point.
pub fn evaluate(expr: *anyopaque, vars: Vector(f64)) f64 {
    return exprtk.evaluate(expr, vars.data.ptr);
}

/// Deinitialize the compiled expression.
pub fn deinit(expr: *anyopaque) void {
    exprtk.deinit(expr);
}

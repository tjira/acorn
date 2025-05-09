//! Wrapper for the tinyexpr C library.

const std = @import("std");

const tinyexpr = @cImport(@cInclude("tinyexpr.h"));

const Vector = @import("vector.zig").Vector;

pub const Expression =     tinyexpr.te_expr;
pub const Variable   = tinyexpr.te_variable;
pub const evaluate   =     tinyexpr.te_eval;

pub fn compile(string: []const u8, vars: []const tinyexpr.te_variable) !*tinyexpr.te_expr {
    var err: i32 = 0;

    const expr = tinyexpr.te_compile(string.ptr, vars.ptr, @intCast(vars.len), &err);

    return if (err == 0) expr else error.InvalidExpression;
}

pub fn makevars(r: Vector(f64), names: []const []const u8, allocator: std.mem.Allocator) ![]tinyexpr.te_variable {
    var vars = try allocator.alloc(tinyexpr.te_variable, names.len);

    for (names, 0..) |name, i| {
        vars[i] = tinyexpr.te_variable{.name = name.ptr, .address = &r.data[i]};
    }

    return vars;
}

pub fn free(expr: *tinyexpr.te_expr, vars: []tinyexpr.te_variable, allocator: std.mem.Allocator) void {
    tinyexpr.te_free(expr); allocator.free(vars);
}

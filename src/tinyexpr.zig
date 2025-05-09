//! Wrapper for the tinyexpr C library.

const tinyexpr = @cImport(@cInclude("tinyexpr.h"));

const std = @import("std");

const Vector = @import("vector.zig").Vector;

pub const Expression =     tinyexpr.te_expr;
pub const Variable   = tinyexpr.te_variable;
pub const evaluate   =     tinyexpr.te_eval;

/// Function to compile a string expression with variables.
pub fn compile(string: []const u8, vars: []const tinyexpr.te_variable) !*tinyexpr.te_expr {
    var err: i32 = 0;

    const expr = tinyexpr.te_compile(string.ptr, vars.ptr, @intCast(vars.len), &err);

    return if (err == 0) expr else error.InvalidExpression;
}

/// Makes the array of variables the tinyexpr library expects.
pub fn makevars(r: Vector(f64), names: []const []const u8, allocator: std.mem.Allocator) ![]tinyexpr.te_variable {
    var vars = try allocator.alloc(tinyexpr.te_variable, names.len + 1);

    const sgn = struct {
        fn fx(x: f64) f64 {return if (x < 0) -1 else 1;}
    };

    vars[0] = tinyexpr.te_variable{.name = "sgn", .address = sgn.fx, .type = tinyexpr.TE_FUNCTION1};

    for (names, 0..) |name, i| {
        vars[i + 1] = tinyexpr.te_variable{.name = name.ptr, .address = &r.data[i]};
    }

    return vars;
}

/// Frees the memory allocated by the expression.
pub fn free(expr: *tinyexpr.te_expr) void {
    tinyexpr.te_free(expr);
}

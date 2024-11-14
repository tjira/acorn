const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn PotentialType(comptime T: type) type {
    return fn (comptime T: type, U: *Matrix(T), r: Vector(T)) void;
}

pub fn doubleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) void {
    U.ptr(0, 0).* = 0.01 * std.math.tanh(0.6 * r.at(0));
    U.ptr(0, 1).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(1, 0).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(1, 1).* = -0.01 * std.math.tanh(0.6 * r.at(0));
}

pub fn evaluate(comptime T: type, U: *Matrix(T), F: *Vector(T), potential: PotentialType(T), r: Vector(T), s: u32, step: T) !void {
    var rp = try r.clone(); var rm = try r.clone(); defer rp.deinit(); defer rm.deinit();
    for (0..r.rows) |i| {
        rp.ptr(i).* += step; potential(T, U, rp); const Up = U.at(s, s);
        rm.ptr(i).* -= step; potential(T, U, rm); const Um = U.at(s, s);
        F.ptr(i).* = -0.5 * (Up - Um) / step;
    }
    potential(T, U, r);
}

pub fn dims(comptime T: type, potential: PotentialType(T)) !u32 {
    if (potential == doubleState1D_1) return 1;
    return error.PotentialNotKnown;
}

pub fn states(comptime T: type, potential: PotentialType(T)) !u32 {
    if (potential == doubleState1D_1) return 2;
    return error.PotentialNotKnown;
}

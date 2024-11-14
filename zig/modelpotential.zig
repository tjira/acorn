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
    // U.ptr(1, 1).* = 0.5 * r.at(0) * r.at(0);
}

pub fn potentialAndForce(comptime T: type, U: *Matrix(T), F: *Vector(T), potential: PotentialType(T), r: Vector(T), s: u8) !void {
    const spatial_step = 0.001;
    for (0..r.rows) |i| {
        const r_plus = try r.clone(); defer r_plus.deinit(); r_plus.ptr(i).* += spatial_step;
        const r_minus = try r.clone(); defer r_minus.deinit(); r_minus.ptr(i).* -= spatial_step;
        var U_plus = try U.clone(); defer U_plus.deinit(); potential(T, &U_plus, r_plus);
        var U_minus = try U.clone(); defer U_minus.deinit(); potential(T, &U_minus, r_minus);
        F.ptr(i).* = -0.5 * (U_plus.at(s, s) - U_minus.at(s, s)) / spatial_step;
    }
    potential(T, U, r);
}

const std = @import("std");

const mat = @import("matrix.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn PotentialType(comptime T: type) type {
    return fn (comptime T: type, U: *Matrix(T), r: Vector(T)) void;
}

pub fn dims(comptime T: type, potential: PotentialType(T)) !u32 {
    if (potential == doubleState1D_1) return 1;
    if (potential == tripleState1D_1) return 1;
    return error.PotentialNotKnown;
}

pub fn states(comptime T: type, potential: PotentialType(T)) !u32 {
    if (potential == doubleState1D_1) return 2;
    if (potential == tripleState1D_1) return 3;
    return error.PotentialNotKnown;
}

pub fn doubleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) void {
    U.ptr(0, 0).* = 0.01 * std.math.tanh(0.6 * r.at(0));
    U.ptr(0, 1).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.01 * std.math.tanh(0.6 * r.at(0));
}

pub fn tripleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) void {
    U.ptr(0, 0).* = 0.03 * (std.math.tanh(1.6 * r.at(0)) + std.math.tanh(1.6 * (r.at(0) + 7)));
    U.ptr(0, 1).* = 0.005 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(0, 2).* = 0.005 * std.math.exp(-(r.at(0) + 7) * (r.at(0) + 7));
    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.03 * (std.math.tanh(1.6 * r.at(0)) + std.math.tanh(1.6 * (r.at(0) - 7)));
    U.ptr(1, 2).* = 0.005 * std.math.exp(-(r.at(0) - 7) * (r.at(0) - 7));
    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = -0.03 * (std.math.tanh(1.6 * (r.at(0) + 7)) - std.math.tanh(1.6 * (r.at(0) - 7)));
}

pub fn evaluate(comptime T: type, potential: PotentialType(T), rvec: Matrix(T), adiabatic: bool) !Matrix(T) {
    var V = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer V.deinit();
    var A = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer A.deinit();
    var C = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer C.deinit();

    var T1 = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer T2.deinit();

    const U = try Matrix(T).init(rvec.rows, V.rows * V.rows, rvec.allocator); var r = try Vector(T).init(rvec.cols, rvec.allocator); defer r.deinit();

    for (0..rvec.rows) |i| {
        
        for (0..rvec.cols) |j| {r.ptr(j).* = rvec.at(i, j);} potential(T, &V, r);

        if (adiabatic) {mat.eigh(T, &A, &C, V, 1e-12, &T1, &T2); @memcpy(V.data, A.data);}

        for (V.data, 0..) |e, j| U.ptr(i, j).* = e;
    }

    return U;
}

pub fn grid(comptime T: type, start: T, end: T, points: u32, dim: u32, allocator: std.mem.Allocator) !Matrix(T) {
    const r = try Matrix(T).init(std.math.pow(u32, points, dim), dim, allocator); const dr = (end - start) / @as(T, @floatFromInt(points - 1));

    for (0..r.rows) |i| for (0..r.cols) |j| {
        r.ptr(i, j).* = start + @as(T, @floatFromInt((i / std.math.pow(usize, points, r.cols - j - 1)) % points)) * dr;
    };

    return r;
}

pub fn write(comptime T: type, path: []const u8, potential: PotentialType(T), start: T, end: T, points: u32, adiabatic: bool, allocator: std.mem.Allocator) !void {
    const r = try grid(T, start, end, points, try dims(T, potential), allocator); defer r.deinit(); const U = try evaluate(T, potential, r, adiabatic); defer U.deinit();

    var V = try Matrix(f64).init(U.rows, r.cols + U.cols, allocator); defer V.deinit(); r.hjoin(&V, U); try V.write(path);
}

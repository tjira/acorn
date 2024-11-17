const std = @import("std");

const mat = @import("matrix.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn PotentialType(comptime T: type) type {
    return fn (comptime T: type, U: *Matrix(T), r: Vector(T)) void;
}

pub fn dims(comptime T: type, potential: PotentialType(T)) !u32 {
    if (potential == harmonic1D_1   ) return 1;
    if (potential == doubleState1D_1) return 1;
    if (potential == tripleState1D_1) return 1;
    return error.PotentialNotKnown;
}

pub fn states(comptime T: type, potential: PotentialType(T)) !u32 {
    if (potential == harmonic1D_1   ) return 1;
    if (potential == doubleState1D_1) return 2;
    if (potential == tripleState1D_1) return 3;
    return error.PotentialNotKnown;
}

pub fn harmonic1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) void {
    U.ptr(0, 0).* = 0.5 * r.at(0) * r.at(0);
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

pub fn evaluate(comptime T: type, U: *Matrix(T), potential: PotentialType(T), rvec: Matrix(T), adiabatic: bool) !void {
    var V = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer V.deinit();
    var A = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer A.deinit();
    var C = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer C.deinit();

    var T1 = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(try states(T, potential), try states(T, potential), rvec.allocator); defer T2.deinit();

    var r = try Vector(T).init(rvec.cols, rvec.allocator); defer r.deinit();

    for (0..rvec.rows) |i| {
        
        for (0..rvec.cols) |j| {r.ptr(j).* = rvec.at(i, j);} potential(T, &V, r);

        if (adiabatic) {mat.eigh(T, &A, &C, V, 1e-12, &T1, &T2); @memcpy(V.data, A.data);}

        for (V.data, 0..) |e, j| U.ptr(i, j).* = e;
    }
}

pub fn kgrid(comptime T: type, k: *Matrix(T), start: T, end: T, points: u32) void {
    k.fill(2 * std.math.pi / asfloat(T, points) / (end - start) * asfloat(T, points - 1));

    for (0..k.rows / points) |i| for (0..points) |j| {k.ptr(i * points + j, k.cols - 1).* *= asfloat(T, j) - asfloat(T, if (j < points / 2) 0 else points);};

    for (0..k.rows) |i| for (0..k.cols) |j| {k.ptr(i, j).* = k.at(i / std.math.pow(usize, points, k.cols - j - 1), k.cols - 1);};
}

pub fn rgrid(comptime T: type, r: *Matrix(T), start: T, end: T, points: u32) void {
    for (0..r.rows) |i| for (0..r.cols) |j| {r.ptr(i, j).* = start + asfloat(T, i / std.math.pow(usize, points, r.cols - j - 1) % points) * (end - start) / asfloat(T, points - 1);};
}

pub fn write(comptime T: type, path: []const u8, potential: PotentialType(T), start: T, end: T, points: u32, adiabatic: bool, allocator: std.mem.Allocator) !void {
    const ndim = try dims(T, potential); const nstate = try states(T, potential); const nrow = std.math.pow(u32, points, ndim); const ncol = nstate * nstate;

    var r = try Matrix(T).init(nrow, ndim, allocator); defer r.deinit(); var U = try Matrix(T).init(nrow, ncol, allocator); defer U.deinit();

    rgrid(T, &r, start, end, points); try evaluate(T, &U, potential, r, adiabatic);

    var V = try Matrix(f64).init(U.rows, r.cols + U.cols, allocator); r.hjoin(&V, U); try V.write(path); V.deinit();
}

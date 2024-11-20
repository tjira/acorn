const std = @import("std"); const Complex = std.math.Complex; const gsl = @cImport(@cInclude("gsl/gsl_eigen.h"));

const mat = @import("matrix.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn PotentialType(comptime T: type) type {
    return fn (comptime T: type, U: *Matrix(T), r: Vector(T)) void;
}

pub fn dims(potential: []const u8) u32 {
    if (std.mem.eql(u8, potential,    "harmonic1D_1")) return 1;
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return 1;
    if (std.mem.eql(u8, potential,       "tully1D_1")) return 1;
    if (std.mem.eql(u8, potential,       "tully1D_2")) return 1;
    return 0;
}

pub fn eval(comptime T: type, U: *Matrix(T), potential: []const u8, r: Vector(T)) void {
    if (std.mem.eql(u8, potential,    "harmonic1D_1"))    harmonic1D_1(T, U, r);
    if (std.mem.eql(u8, potential, "doubleState1D_1")) doubleState1D_1(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_1")) tripleState1D_1(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_1"))       tully1D_1(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_2"))       tully1D_2(T, U, r);
}

pub fn states(potential: []const u8) u32 {
    if (std.mem.eql(u8, potential,    "harmonic1D_1")) return 1;
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return 2;
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return 3;
    if (std.mem.eql(u8, potential,       "tully1D_1")) return 2;
    if (std.mem.eql(u8, potential,       "tully1D_2")) return 2;
    return 0;
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

pub fn tully1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) void {
    U.ptr(0, 0).* = if (r.at(0) > 0) 0.01 * (1 - std.math.exp(-1.6 * r.at(0))) else -0.01 * (1 - std.math.exp(1.6 * r.at(0)));
    U.ptr(0, 1).* = 0.005 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = if (r.at(0) > 0) -0.01 * (1 - std.math.exp(-1.6 * r.at(0))) else 0.01 * (1 - std.math.exp(1.6 * r.at(0)));
}

pub fn tully1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) void {
    U.ptr(0, 0).* = 0;
    U.ptr(0, 1).* = 0.0150 * std.math.exp(-0.06 * r.at(0) * r.at(0));
    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.1 * std.math.exp(-0.28 * r.at(0) * r.at(0)) + 0.05;
}

pub fn kgrid(comptime T: type, k: *Matrix(T), start: T, end: T, points: u32) void {
    k.fill(2 * std.math.pi / asfloat(T, points) / (end - start) * asfloat(T, points - 1));

    for (0..k.rows / points) |i| for (0..points) |j| {
        k.ptr(i * points + j, k.cols - 1).* *= asfloat(T, j) - asfloat(T, if (j < points / 2) 0 else points);
    };

    for (0..k.rows) |i| for (0..k.cols) |j| {
        k.ptr(i, j).* = k.at(i / std.math.pow(usize, points, k.cols - j - 1), k.cols - 1);
    };
}

pub fn rgrid(comptime T: type, r: *Matrix(T), start: T, end: T, points: u32) void {
    for (0..r.rows) |i| for (0..r.cols) |j| {
        r.ptr(i, j).* = start + asfloat(T, i / std.math.pow(usize, points, r.cols - j - 1) % points) * (end - start) / asfloat(T, points - 1);
    };
}

pub fn write(comptime T: type, path: []const u8, potential: []const u8, start: T, end: T, points: u32, adiabatic: bool, allocator: std.mem.Allocator) !void {
    const GSLW = gsl.gsl_eigen_symmv_alloc(states(potential)); defer gsl.gsl_eigen_symmv_free(GSLW);

    var U  = try Matrix(T).init(states(potential), states(potential), allocator); defer  U.deinit();
    var UA = try Matrix(T).init(states(potential), states(potential), allocator); defer UA.deinit();
    var UC = try Matrix(T).init(states(potential), states(potential), allocator); defer UC.deinit();

    var R = try Matrix(T).init(std.math.pow(u32, points, dims(potential)), dims(potential)                      , allocator); defer R.deinit();
    var V = try Matrix(T).init(std.math.pow(u32, points, dims(potential)), states(potential) * states(potential), allocator); defer V.deinit();

    rgrid(T, &R, start, end, points); UA.fill(0);

    for (0..R.rows) |i| {eval(T, &U, potential, R.rowptr(i).vectorptr()); if (adiabatic) {mat.eigh(T, &UA, &UC, U, GSLW); @memcpy(U.data, UA.data);} for (U.data, 0..) |e, j| V.ptr(i, j).* = e;}

    var VT = try Matrix(T).init(V.rows, R.cols + V.cols, allocator); mat.hjoin(T, &VT, R, V); try VT.write(path); VT.deinit();
}

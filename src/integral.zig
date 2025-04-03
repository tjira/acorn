//! Integral engine.

const std = @import("std");

const ContractedGaussian  = @import("contractedgaussian.zig" ).ContractedGaussian;
const Matrix              = @import("matrix.zig"             ).Matrix            ;
const Tensor              = @import("tensor.zig"             ).Tensor            ;
const System              = @import("system.zig"             ).System            ;

/// Compute the kinetic matrix.
pub fn coulomb(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), allocator: std.mem.Allocator) !Tensor(T) {
    var J = try Tensor(T).init(&[_]usize{basis.items.len, basis.items.len, basis.items.len, basis.items.len}, allocator);

    for (0..J.shape[0]) |i| for (i..J.shape[1]) |j| for (i..J.shape[2]) |k| for ((if (i == k) j else k)..J.shape[3]) |l| {

        const I = basis.items[i].coulomb(basis.items[j], basis.items[k], basis.items[l]);

        J.ptr(&[_]usize{i, j, k, l}).* = I;
        J.ptr(&[_]usize{i, j, l, k}).* = I;
        J.ptr(&[_]usize{j, i, k, l}).* = I;
        J.ptr(&[_]usize{j, i, l, k}).* = I;

        J.ptr(&[_]usize{k, l, i, j}).* = I;
        J.ptr(&[_]usize{k, l, j, i}).* = I;
        J.ptr(&[_]usize{l, k, i, j}).* = I;
        J.ptr(&[_]usize{l, k, j, i}).* = I;
    };

    return J;
}

/// Compute the kinetic matrix.
pub fn kinetic(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), allocator: std.mem.Allocator) !Matrix(T) {
    var K = try Matrix(T).init(basis.items.len, basis.items.len, allocator);

    for (0..K.rows) |i| for (i..K.cols) |j| {
        K.ptr(i, j).* = basis.items[i].kinetic(basis.items[j]); K.ptr(j, i).* = K.at(i, j);
    };

    return K;
}

/// Compute the kinetic matrix.
pub fn nuclear(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), system: System(T), allocator: std.mem.Allocator) !Matrix(T) {
    var V = try Matrix(T).init(basis.items.len, basis.items.len, allocator);

    for (0.. V.rows) |i| for (i..V.cols) |j| {
        V.ptr(i, j).* = basis.items[i].nuclear(basis.items[j], system); V.ptr(j, i).* = V.at(i, j);
    };

    return V;
}

/// Compute the overlap matrix.
pub fn overlap(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), allocator: std.mem.Allocator) !Matrix(T) {
    var S = try Matrix(T).init(basis.items.len, basis.items.len, allocator); S.fill(0);

    for (0..S.rows) |i| for (i..S.cols) |j| {
        S.ptr(i, j).* = basis.items[i].overlap(basis.items[j]); S.ptr(j, i).* = S.at(i, j);
    };

    return S;
}

//! Integral engine.

const std = @import("std");

const ContractedGaussian  = @import("contractedgaussian.zig" ).ContractedGaussian;
const Matrix              = @import("matrix.zig"             ).Matrix            ;
const Tensor              = @import("tensor.zig"             ).Tensor            ;
const System              = @import("system.zig"             ).System            ;

/// Compute the kinetic matrix.
pub fn coulomb(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), allocator: std.mem.Allocator) !Tensor(T) {
    var J = try Tensor(T).init(&[_]usize{basis.items.len, basis.items.len, basis.items.len, basis.items.len}, allocator);

    for (basis.items, 0..) |cg1, i| for (basis.items, 0..) |cg2, j| for (basis.items, 0..) |cg3, k| for (basis.items, 0..) |cg4, l| {
        J.ptr(&[_]usize{i, j, k, l}).* += cg1.coulomb(cg2, cg3, cg4);
    };

    return J;
}

/// Compute the kinetic matrix.
pub fn kinetic(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), allocator: std.mem.Allocator) !Matrix(T) {
    var K = try Matrix(f64).init(basis.items.len, basis.items.len, allocator);

    for (basis.items, 0..) |cg1, i| for (basis.items, 0..) |cg2, j| {
        K.ptr(i, j).* += cg1.kinetic(cg2);
    };

    return K;
}

/// Compute the kinetic matrix.
pub fn nuclear(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), system: System(T), allocator: std.mem.Allocator) !Matrix(T) {
    var V = try Matrix(T).init(basis.items.len, basis.items.len, allocator);

    for (basis.items, 0..) |cg1, i| for (basis.items, 0..) |cg2, j| {
        V.ptr(i, j).* += cg1.nuclear(cg2, system);
    };

    return V;
}

/// Compute the overlap matrix.
pub fn overlap(comptime T: type, basis: std.ArrayList(ContractedGaussian(T)), allocator: std.mem.Allocator) !Matrix(T) {
    var S = try Matrix(T).init(basis.items.len, basis.items.len, allocator); S.fill(0);

    for (basis.items, 0..) |cg1, i| for (basis.items, 0..) |cg2, j| {
        S.ptr(i, j).* += cg1.overlap(cg2);
    };

    return S;
}

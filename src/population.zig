//! File that contains functions for population analysis.

const std = @import("std");

const mat = @import("matrix.zig");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

/// Mulliken population analysis. The function outputs the Mulliken population of each atom in the system given the system struct, overlap matrix and density matrix.
pub fn mulliken(comptime T: type, system: System(T), basis: std.ArrayList(T), S_AO: Matrix(T), D_AO: Matrix(T), allocator: std.mem.Allocator) !Vector(T) {
    var m = try Vector(T).init(system.atoms.rows, allocator);

    for (0..system.atoms.rows) |i| m.ptr(i).* = system.atoms.at(i);

    const bf2atom = try system.getBasisMap(basis, allocator); defer bf2atom.deinit();

    for (0..S_AO.rows) |i| for (0..S_AO.cols) |j| {
        m.ptr(bf2atom.items[i % bf2atom.items.len]).* -= D_AO.at(j, i) * S_AO.at(j, i);
    };

    return m;
}

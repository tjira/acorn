const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn ModelPotential(comptime T: type) type {
    return struct {
        pub fn dim(potential: fn (r: Vector(T)) anyerror!Matrix(T)) !usize {
            if (potential == ModelPotential(T).tully1    ) return 1;
            if (potential == ModelPotential(T).harmonic1D) return 1;
            return error.PotentialNotKnown;
        }
        pub fn tully1(r: Vector(T)) !Matrix(T) {
            const U = try Matrix(T).zero(2, 2, r.data.allocator);
            if (r.at(0) >= 0) U.set(0, 0,  0.01 * (1 - std.math.exp(-1.6 * r.at(0))));
            if (r.at(0) <  0) U.set(0, 0, -0.01 * (1 - std.math.exp( 1.6 * r.at(0))));
            U.set(0, 1, 0.005 * std.math.exp(-r.at(0) * r.at(0)));
            U.set(1, 0, 0.005 * std.math.exp(-r.at(0) * r.at(0)));
            if (r.at(0) >= 0) U.set(1, 1, -0.01 * (1 - std.math.exp(-1.6 * r.at(0))));
            if (r.at(0) <  0) U.set(1, 1,  0.01 * (1 - std.math.exp( 1.6 * r.at(0))));
            return U;
        }
        pub fn harmonic1D(r: Vector(T)) !Matrix(T) {
            const U = try Matrix(T).zero(1, 1, r.data.allocator);
            try U.set(0, 0, 0.5 * r.at(0) * r.at(0));
            return U;
        }
    };
}

pub fn evaluateModelPotential(comptime T: type, potential: fn (r: Vector(T)) anyerror!Matrix(T), from: T, to: T, points: usize, allocator: std.mem.Allocator) !Matrix(T) {
    const dim = try ModelPotential(T).dim(potential);
    const result = try Matrix(T).zero(std.math.pow(usize, points, dim), dim + (try potential(try Vector(T).zero(dim, allocator))).data.items.len, allocator);
    const dr = (to - from) / @as(T, @floatFromInt(points - 1));
    for (0..result.rows) |i| {
        var r = try Vector(T).zero(dim, allocator);
        for (0..dim) |j| {
            r.set(j, from + @as(T, @floatFromInt((i / std.math.pow(usize, points, dim - j - 1)) % points)) * dr);
            result.set(i, j, r.at(j));
        }
        const U = try potential(r);
        for (0..U.data.items.len) |j| result.set(i, dim + j, U.data.items[j]);
    }
    return result;
}

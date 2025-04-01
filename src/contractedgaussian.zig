//! Gaussian type definitions.

const std = @import("std");

const System = @import("system.zig").System;

const sum = @import("helper.zig").sum;

/// Contracted Gaussian type.
pub fn ContractedGaussian(comptime T: type) type {
    return struct {
        A: [3]T, a: [3]T, c: []T, alpha: []T, l: u8,

        /// Initialize a contracted Gaussian.
        pub fn init(A: [3]T, a: [3]T, c: []const T, alpha: []const T, allocator: std.mem.Allocator) !ContractedGaussian(T) {

            const g = ContractedGaussian(T) {.A = A, .a = a, .c = try allocator.alloc(T, c.len), .alpha = try allocator.alloc(T, alpha.len), .l = @as(u8, @intFromFloat(sum(T, &a)))};

            @memcpy(g.c, c); @memcpy(g.alpha, alpha);

            for (g.c, 0..) |*ci, i| ci.* *= std.math.pow(T, 2 * g.alpha[i] / std.math.pi, 0.75);

            return g;
        }

        /// Deinitialize a contracted Gaussian.
        pub fn deinit(self: ContractedGaussian(T), allocator: std.mem.Allocator) void {
            allocator.free(self.c); allocator.free(self.alpha);
        }

        /// Compute the Coulomb integral between four contracted Gaussians.
        pub fn coulomb(self: ContractedGaussian(T), other1: ContractedGaussian(T), other2: ContractedGaussian(T), other3: ContractedGaussian(T)) !T {
            var e: T = 0;

            for (self.c, 0..) |ci, i| for (other1.c, 0..) |cj, j| for (other2.c, 0..) |ck, k| for (other3.c, 0..) |cl, l| {

                if (self.l == 0 and other1.l == 0 and other2.l == 0 and other3.l == 0) {e += ci * cj * ck * cl * self.coulombPrimitiveS(other1, other2, other3, i, j, k, l);}

                else return error.UnsupportedAngularMomentum;
            };

            return e;
        }

        /// Compute the Coulomb integral between four primitive Gaussians.
        pub fn coulombPrimitiveS(self: ContractedGaussian(T), other1: ContractedGaussian(T), other2: ContractedGaussian(T), other3: ContractedGaussian(T), i: usize, j: usize, k: usize, l: usize) T {
            _ = self; _ = other1; _ = other2; _ = other3; _ = i; _ = j; _ = k; _ = l; return 0;
        }

        /// Compute the kinetic integral between two contracted Gaussians.
        pub fn kinetic(self: ContractedGaussian(T), other: ContractedGaussian(T)) !T {
            var t: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                if (self.l == 0 and other.l == 0) {t += ci * cj * self.kineticPrimitiveS(other, i, j);}

                else return error.UnsupportedAngularMomentum;
            };

            return t;
        }

        /// Compute the kinetic integral between two primitive Gaussians.
        pub fn kineticPrimitiveS(self: ContractedGaussian(T), other: ContractedGaussian(T), i: usize, j: usize) T {
            _ = self; _ = other; _ = i; _ = j; return 0;
        }

        /// Compute the nuclear integral between two contracted Gaussians.
        pub fn nuclear(self: ContractedGaussian(T), other: ContractedGaussian(T), system: System(T)) !T {
            var v: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                if (self.l == 0 and other.l == 0) {v += ci * cj * self.nuclearPrimitiveS(other, i, j, system);}

                else return error.UnsupportedAngularMomentum;
            };

            return v;
        }

        /// Compute the nuclear integral between two primitive Gaussians.
        pub fn nuclearPrimitiveS(self: ContractedGaussian(T), other: ContractedGaussian(T), i: usize, j: usize, system: System(T)) T {
            _ = self; _ = other; _ = i; _ = j; _ = system; return 0;
        }

        /// Compute the overlap integral between two contracted Gaussians.
        pub fn overlap(self: ContractedGaussian(T), other: ContractedGaussian(T)) !T {
            var s: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                if (self.l == 0 and other.l == 0) {s += ci * cj * self.overlapPrimitiveS(other, i, j);}

                else return error.UnsupportedAngularMomentum;
            };

            return s;
        }

        /// Compute the overlap integral between two primitive Gaussians.
        pub fn overlapPrimitiveS(self: ContractedGaussian(T), other: ContractedGaussian(T), i: usize, j: usize) T {
            const factor = std.math.pow(T, std.math.pi / (self.alpha[i] + other.alpha[j]), 1.5);

            const vectorsq = std.math.pow(T, self.A[0] - other.A[0], 2) + std.math.pow(T, self.A[1] - other.A[1], 2) + std.math.pow(T, self.A[2] - other.A[2], 2);

            const exponent = -self.alpha[i] * other.alpha[j] * vectorsq / (self.alpha[i] + other.alpha[j]);

            return factor * std.math.exp(exponent);
        }
    };
}

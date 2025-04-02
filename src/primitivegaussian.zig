//! Primitive Gaussian type definitions.

const std = @import("std");

const System = @import("system.zig").System;

const sum = @import("helper.zig").sum;

/// Primitive Gaussian type.
pub fn PrimitiveGaussian(comptime T: type) type {
    return struct {
        A: [3]T, a: [3]T, alpha: T, l: i8,

        /// Compute the coulomb integral between four primitive Gaussians.
        pub fn coulomb(self: PrimitiveGaussian(T), other1: PrimitiveGaussian(T), other2: PrimitiveGaussian(T), other3: PrimitiveGaussian(T)) !T {
            _ = self; _ = other1; _ = other2; _ = other3; return 0;
        }

        /// Compute the kinetic integral between two primitive Gaussians.
        pub fn kinetic(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) !T {
            _ = self; _ = other; return 0;
        }

        /// Compute the nuclear integral between two primitive Gaussians.
        pub fn nuclear(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T), system: System(T)) !T {
            _ = self; _ = other; _ = system; return 0;
        }

        /// Compute the overlap integral between two primitive Gaussians.
        pub fn overlap(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) !T {
            const p = self.alpha + other.alpha; const mu = self.alpha * other.alpha / p;

            const XAB: [3]T = .{
                self.A[0] - other.A[0],
                self.A[1] - other.A[1],
                self.A[2] - other.A[2]
            };

            const XPA: [3]T = .{
                (self.alpha * self.A[0] + other.alpha * other.A[0]) / p - self.A[0],
                (self.alpha * self.A[1] + other.alpha * other.A[1]) / p - self.A[1],
                (self.alpha * self.A[2] + other.alpha * other.A[2]) / p - self.A[2]
            };

            const XPB: [3]T = .{
                (self.alpha * self.A[0] + other.alpha * other.A[0]) / p - other.A[0],
                (self.alpha * self.A[1] + other.alpha * other.A[1]) / p - other.A[1],
                (self.alpha * self.A[2] + other.alpha * other.A[2]) / p - other.A[2]
            };

            const Sij = try overlapSingle(.{self.a[0], other.a[0]}, p, mu, XAB[0], XPA[0], XPB[0]);
            const Skl = try overlapSingle(.{self.a[1], other.a[1]}, p, mu, XAB[1], XPA[1], XPB[1]);
            const Smn = try overlapSingle(.{self.a[2], other.a[2]}, p, mu, XAB[2], XPA[2], XPB[2]);

            return Sij * Skl * Smn;
        }

        /// Compute the overlap integral between two primitive Gaussians in one dimension using the Obara-Saika recursion.
        pub fn overlapSingle(ij: [2]T, p: T, mu: T, XAB: T, XPA: T, XPB: T) !T {
            if (ij[0] == 0 and ij[1] == 0) {return std.math.sqrt(std.math.pi / p) * std.math.exp(-mu * XAB * XAB);}

            else if (ij[0] > 0) {

                const Sij0 = try overlapSingle(.{ij[0] - 1, ij[1]    }, p, mu, XAB, XPA, XPB);
                const Sij1 = try overlapSingle(.{ij[0] - 2, ij[1]    }, p, mu, XAB, XPA, XPB);
                const Sij2 = try overlapSingle(.{ij[0] - 1, ij[1] - 1}, p, mu, XAB, XPA, XPB);

                return XPA * Sij0 + 1 / (2 * p) * ((ij[0] - 1) * Sij1 + ij[1] * Sij2);
            }

            else if (ij[1] > 0) {

                const Sij0 = try overlapSingle(.{ij[0]    , ij[1] - 1}, p, mu, XAB, XPA, XPB);
                const Sij1 = try overlapSingle(.{ij[0] - 1, ij[1] - 1}, p, mu, XAB, XPA, XPB);
                const Sij2 = try overlapSingle(.{ij[0]    , ij[1] - 2}, p, mu, XAB, XPA, XPB);

                return XPB * Sij0 + 1 / (2 * p) * (ij[0] * Sij1 + (ij[1] - 1) * Sij2);
            }

            else return 0;
        }
    };
}

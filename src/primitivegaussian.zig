//! Primitive Gaussian type definitions.

const std = @import("std");

const System = @import("system.zig").System;

const sum = @import("helper.zig").sum;

/// Primitive Gaussian type.
pub fn PrimitiveGaussian(comptime T: type) type {
    return struct {
        A: [3]T, a: [3]T, alpha: T, l: i8,

        /// Compute the coulomb integral between four primitive Gaussians using the Obara-Saika recursion.
        pub fn coulomb(self: PrimitiveGaussian(T), other1: PrimitiveGaussian(T), other2: PrimitiveGaussian(T), other3: PrimitiveGaussian(T)) !T {
            _ = self; _ = other1; _ = other2; _ = other3; return 0;
        }

        /// Compute the kinetic integral between two primitive Gaussians using the Obara-Saika recursion.
        pub fn kinetic(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) !T {
            _ = self; _ = other; return 0;
        }

        /// Compute the nuclear integral between two primitive Gaussians using the Obara-Saika recursion.
        pub fn nuclear(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T), system: System(T)) !T {
            _ = self; _ = other; _ = system; return 0;
        }

        /// Compute the overlap integral between two primitive Gaussians using the Obara-Saika recursion.
        pub fn overlap(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) !T {
            var Sij: T = 0; var Skl: T = 0; var Smn: T = 0; const p = self.alpha + other.alpha; const mu = self.alpha * other.alpha / p;

            const XAB: [3]T = .{
                other.A[0] - self.A[0],
                other.A[1] - self.A[1],
                other.A[2] - self.A[2]
            };

            const XPA: [3]T = .{
                self.A[0] - (self.alpha * self.A[0] + other.alpha * other.A[0]) / p,
                self.A[1] - (self.alpha * self.A[1] + other.alpha * other.A[1]) / p,
                self.A[2] - (self.alpha * self.A[2] + other.alpha * other.A[2]) / p
            };

            const XPB: [3]T = .{
                other.A[0] - (self.alpha * self.A[0] + other.alpha * other.A[0]) / p,
                other.A[1] - (self.alpha * self.A[1] + other.alpha * other.A[1]) / p,
                other.A[2] - (self.alpha * self.A[2] + other.alpha * other.A[2]) / p
            };

            _ = XPA; _ = XPB;

            if (self.l == 0 and other.l == 0) {
                Sij = std.math.sqrt(std.math.pi / p) * std.math.exp(-mu * XAB[0] * XAB[0]);
                Skl = std.math.sqrt(std.math.pi / p) * std.math.exp(-mu * XAB[1] * XAB[1]);
                Smn = std.math.sqrt(std.math.pi / p) * std.math.exp(-mu * XAB[2] * XAB[2]);
            }

            else return error.InvalidAngularMomentum; return Sij * Skl * Smn;
        }
    };
}

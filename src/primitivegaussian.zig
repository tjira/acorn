//! Primitive Gaussian type definitions.

const std = @import("std");

const System = @import("system.zig").System;

const dfact = @import("helper.zig").dfact;
const powi  = @import("helper.zig").powi ;
const sum   = @import("helper.zig").sum  ;

/// Primitive Gaussian type.
pub fn PrimitiveGaussian(comptime T: type) type {
    return struct {
        A: [3]T, a: [3]T, alpha: T, l: i8,

        /// Compute the coulomb integral between four primitive Gaussians.
        pub fn coulomb(self: PrimitiveGaussian(T), other1: PrimitiveGaussian(T), other2: PrimitiveGaussian(T), other3: PrimitiveGaussian(T)) T {
            _ = self; _ = other1; _ = other2; _ = other3; return 0;
        }

        /// Compute the kinetic integral between two primitive Gaussians.
        pub fn kinetic(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) T {
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

            const Sij = overlapSingle(.{self.a[0], other.a[0]}, p, mu, XAB[0], XPA[0], XPB[0]);
            const Skl = overlapSingle(.{self.a[1], other.a[1]}, p, mu, XAB[1], XPA[1], XPB[1]);
            const Smn = overlapSingle(.{self.a[2], other.a[2]}, p, mu, XAB[2], XPA[2], XPB[2]);

            const Tij = kineticSingle(.{self.a[0], other.a[0]}, .{self.alpha, other.alpha}, p, mu, XAB[0], XPA[0], XPB[0]);
            const Tkl = kineticSingle(.{self.a[1], other.a[1]}, .{self.alpha, other.alpha}, p, mu, XAB[1], XPA[1], XPB[1]);
            const Tmn = kineticSingle(.{self.a[2], other.a[2]}, .{self.alpha, other.alpha}, p, mu, XAB[2], XPA[2], XPB[2]);

            return Tij * Skl * Smn + Sij * Tkl * Smn + Sij * Skl * Tmn;
        }

        /// Compute the kinetic integral between two primitive Gaussians in one dimension using the Obara-Saika recursion.
        pub fn kineticSingle(ij: [2]T, ab: [2]T, p: T, mu: T, XAB: T, XPA: T, XPB: T) T {
            if (ij[0] == 0 and ij[1] == 0) {
                return (ab[0] - 2 * ab[0] * ab[0] * (XPA * XPA + 1 / (2 * p))) * overlapSingle(.{ij[0], ij[1]}, p, mu, XAB, XPA, XPB);
            }

            if (ij[0] > 0) {

                const Tij0 = kineticSingle(.{ij[0] - 1, ij[1]    }, ab, p, mu, XAB, XPA, XPB);
                const Tij1 = kineticSingle(.{ij[0] - 2, ij[1]    }, ab, p, mu, XAB, XPA, XPB);
                const Tij2 = kineticSingle(.{ij[0] - 1, ij[1] - 1}, ab, p, mu, XAB, XPA, XPB);

                const Sij0 = overlapSingle(.{ij[0]    , ij[1]    }, p, mu, XAB, XPA, XPB);
                const Sij1 = overlapSingle(.{ij[0] - 2, ij[1]    }, p, mu, XAB, XPA, XPB);

                return XPA * Tij0 + 1 / (2 * p) * ((ij[0] - 1) * Tij1 + ij[1] * Tij2) + ab[1] / p * (2 * ab[0] * Sij0 - (ij[0] - 1) * Sij1);
            }

            if (ij[1] > 0) {

                const Tij0 = kineticSingle(.{ij[0]    , ij[1] - 1}, ab, p, mu, XAB, XPA, XPB);
                const Tij1 = kineticSingle(.{ij[0] - 1, ij[1] - 1}, ab, p, mu, XAB, XPA, XPB);
                const Tij2 = kineticSingle(.{ij[0]    , ij[1] - 2}, ab, p, mu, XAB, XPA, XPB);

                const Sij0 = overlapSingle(.{ij[0]    , ij[1]    }, p, mu, XAB, XPA, XPB);
                const Sij1 = overlapSingle(.{ij[0]    , ij[1] - 2}, p, mu, XAB, XPA, XPB);

                return XPB * Tij0 + 1 / (2 * p) * (ij[0] * Tij1 + (ij[1] - 1) * Tij2) + ab[0] / p * (2 * ab[1] * Sij0 - (ij[1] - 1) * Sij1);
            }

            else return 0;
        }

        /// Calculate the norm of the primitive gaussian.
        pub fn norm(self: PrimitiveGaussian(T)) T {
            const Nij = dfact(2 * self.a[0] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / powi(4 * self.alpha, @as(u32, @intFromFloat(self.a[0])));
            const Nkl = dfact(2 * self.a[1] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / powi(4 * self.alpha, @as(u32, @intFromFloat(self.a[1])));
            const Nmn = dfact(2 * self.a[2] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / powi(4 * self.alpha, @as(u32, @intFromFloat(self.a[2])));

            return std.math.sqrt(Nij * Nkl * Nmn);
        }

        /// Compute the nuclear integral between two primitive Gaussians.
        pub fn nuclear(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T), system: System(T)) T {
            _ = self; _ = other; _ = system; return 0;
        }

        /// Compute the overlap integral between two primitive Gaussians.
        pub fn overlap(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) T {
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

            const Sij = overlapSingle(.{self.a[0], other.a[0]}, p, mu, XAB[0], XPA[0], XPB[0]);
            const Skl = overlapSingle(.{self.a[1], other.a[1]}, p, mu, XAB[1], XPA[1], XPB[1]);
            const Smn = overlapSingle(.{self.a[2], other.a[2]}, p, mu, XAB[2], XPA[2], XPB[2]);

            return Sij * Skl * Smn;
        }

        /// Compute the overlap integral between two primitive Gaussians in one dimension using the Obara-Saika recursion.
        pub fn overlapSingle(ij: [2]T, p: T, mu: T, XAB: T, XPA: T, XPB: T) T {
            if (ij[0] == 0 and ij[1] == 0) {return std.math.sqrt(std.math.pi / p) * std.math.exp(-mu * XAB * XAB);}

            else if (ij[0] > 0) {

                const Sij0 = overlapSingle(.{ij[0] - 1, ij[1]    }, p, mu, XAB, XPA, XPB);
                const Sij1 = overlapSingle(.{ij[0] - 2, ij[1]    }, p, mu, XAB, XPA, XPB);
                const Sij2 = overlapSingle(.{ij[0] - 1, ij[1] - 1}, p, mu, XAB, XPA, XPB);

                return XPA * Sij0 + 1 / (2 * p) * ((ij[0] - 1) * Sij1 + ij[1] * Sij2);
            }

            else if (ij[1] > 0) {

                const Sij0 = overlapSingle(.{ij[0]    , ij[1] - 1}, p, mu, XAB, XPA, XPB);
                const Sij1 = overlapSingle(.{ij[0] - 1, ij[1] - 1}, p, mu, XAB, XPA, XPB);
                const Sij2 = overlapSingle(.{ij[0]    , ij[1] - 2}, p, mu, XAB, XPA, XPB);

                return XPB * Sij0 + 1 / (2 * p) * (ij[0] * Sij1 + (ij[1] - 1) * Sij2);
            }

            else return 0;
        }
    };
}

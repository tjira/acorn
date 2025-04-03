//! Primitive Gaussian type definitions.

const std = @import("std");

const mth = @import("math.zig");

const System = @import("system.zig").System;

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
            const pgpl = PrimitiveGaussian(T){.A = other.A, .a = .{other.a[0] + 2, other.a[1], other.a[2]}, .alpha = other.alpha, .l = other.l + 2};
            const pgpm = PrimitiveGaussian(T){.A = other.A, .a = .{other.a[0], other.a[1] + 2, other.a[2]}, .alpha = other.alpha, .l = other.l + 2};
            const pgpn = PrimitiveGaussian(T){.A = other.A, .a = .{other.a[0], other.a[1], other.a[2] + 2}, .alpha = other.alpha, .l = other.l + 2};
            const pgml = PrimitiveGaussian(T){.A = other.A, .a = .{other.a[0] - 2, other.a[1], other.a[2]}, .alpha = other.alpha, .l = other.l - 2};
            const pgmm = PrimitiveGaussian(T){.A = other.A, .a = .{other.a[0], other.a[1] - 2, other.a[2]}, .alpha = other.alpha, .l = other.l - 2};
            const pgmn = PrimitiveGaussian(T){.A = other.A, .a = .{other.a[0], other.a[1], other.a[2] - 2}, .alpha = other.alpha, .l = other.l - 2};

            const T0 = other.alpha * (2 * mth.sum(T, &other.a) + 3) * self.overlap(other);

            const T1 = -2 * std.math.pow(T, other.alpha, 2.0) * (self.overlap(pgpl) + self.overlap(pgpm) + self.overlap(pgpn));

            const T2 = -0.5 * (other.a[0] * (other.a[0] - 1) * self.overlap(pgml) + other.a[1] * (other.a[1] - 1) * self.overlap(pgmm) + other.a[2] * (other.a[2] - 1) * self.overlap(pgmn));

            return T0 + T1 + T2;
        }

        /// Calculate the norm of the primitive gaussian.
        pub fn norm(self: PrimitiveGaussian(T)) T {
            const Nij = mth.dfact(2 * self.a[0] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / mth.powi(4 * self.alpha, @as(u32, @intFromFloat(self.a[0])));
            const Nkl = mth.dfact(2 * self.a[1] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / mth.powi(4 * self.alpha, @as(u32, @intFromFloat(self.a[1])));
            const Nmn = mth.dfact(2 * self.a[2] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / mth.powi(4 * self.alpha, @as(u32, @intFromFloat(self.a[2])));

            return std.math.sqrt(Nij * Nkl * Nmn);
        }

        /// Compute the nuclear integral between two primitive Gaussians.
        pub fn nuclear(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T), system: System(T)) T {
            const p = self.alpha + other.alpha; var n: T = 0;

            const XAB: [3]T = .{
                self.A[0] - other.A[0],
                self.A[1] - other.A[1],
                self.A[2] - other.A[2]
            };

            for (0..system.atoms.rows) |i| {

                const XPC: [3]T = .{
                    (self.alpha * self.A[0] + other.alpha * other.A[0]) / p - system.coords.at(i, 0),
                    (self.alpha * self.A[1] + other.alpha * other.A[1]) / p - system.coords.at(i, 1),
                    (self.alpha * self.A[2] + other.alpha * other.A[2]) / p - system.coords.at(i, 2)
                };

                for (0..@as(usize, @intFromFloat(self.a[0] + other.a[0])) + 1) |t| {
                    for (0..@as(usize, @intFromFloat(self.a[1] + other.a[1])) + 1) |u| {
                        for (0..@as(usize, @intFromFloat(self.a[2] + other.a[2])) + 1) |v| {

                            const Eij = hermitec(.{self.a[0], other.a[0]}, .{self.alpha, other.alpha}, XAB[0], @floatFromInt(t));
                            const Ekl = hermitec(.{self.a[1], other.a[1]}, .{self.alpha, other.alpha}, XAB[1], @floatFromInt(u));
                            const Emn = hermitec(.{self.a[2], other.a[2]}, .{self.alpha, other.alpha}, XAB[2], @floatFromInt(v));

                            n -= system.atoms.at(i) * Eij * Ekl * Emn * hermiteint([3]T{@floatFromInt(t), @floatFromInt(u), @floatFromInt(v)}, XPC, p, 0);
                        }
                    }
                }
            }

            return 2 * std.math.pi / p * n;
        }

        /// Compute the Hermite Gaussian coefficients in one dimension.
        pub fn hermitec(ij: [2]T, ab: [2]T, Q: T, t: T) T {
            const p = ab[0] + ab[1]; const q = ab[0] * ab[1] / p;

            if (t < 0 or t > ij[0] + ij[1]) {return 0;}

            else if (ij[0] == 0 and ij[1] == 0 and t == 0) {return std.math.exp(-q * Q * Q);}

            else if (ij[1] == 0) {

                const E1 = hermitec(.{ij[0] - 1, ij[1]}, ab, Q, t - 1);
                const E2 = hermitec(.{ij[0] - 1, ij[1]}, ab, Q, t    );
                const E3 = hermitec(.{ij[0] - 1, ij[1]}, ab, Q, t + 1);

                return (1 / (2 * p)) * E1 - (q * Q / ab[0]) * E2 + (t + 1) * E3;
            }

            const E1 = hermitec(.{ij[0], ij[1] - 1}, ab, Q, t - 1);
            const E2 = hermitec(.{ij[0], ij[1] - 1}, ab, Q, t    );
            const E3 = hermitec(.{ij[0], ij[1] - 1}, ab, Q, t + 1);

            return (1 / (2 * p)) * E1 + (q * Q / ab[1]) * E2 + (t + 1) * E3;
        }

        /// Compute the Hermite integrals.
        pub fn hermiteint(tuv: [3]T, XPC: [3]T, p: T, n: T) T {
            var I: T = 0;

            if (tuv[0] == 0 and tuv[1] == 0 and tuv[2] == 0) {
                I += std.math.pow(T, -2 * p, n) * mth.boys(p * (XPC[0] * XPC[0] + XPC[1] * XPC[1] + XPC[2] * XPC[2]), n);
            }

            else if (tuv[0] == 0 and tuv[1] == 0) {
                if (tuv[2] > 1) {
                    I += (tuv[2] - 1) * hermiteint(.{tuv[0], tuv[1], tuv[2] - 2}, XPC, p, n + 1);
                }
                I += XPC[2] * hermiteint(.{tuv[0], tuv[1], tuv[2] - 1}, XPC, p, n + 1);
            }

            else if (tuv[0] == 0) {
                if (tuv[1] > 1) {
                    I += (tuv[1] - 1) * hermiteint(.{tuv[0], tuv[1] - 2, tuv[2]}, XPC, p, n + 1);
                }
                I += XPC[1] * hermiteint(.{tuv[0], tuv[1] - 1, tuv[2]}, XPC, p, n + 1);
            }

            else {
                if (tuv[0] > 1) {
                    I += (tuv[0] - 1) * hermiteint(.{tuv[0] - 2, tuv[1], tuv[2]}, XPC, p, n + 1);
                }
                I += XPC[0] * hermiteint(.{tuv[0] - 1, tuv[1], tuv[2]}, XPC, p, n + 1);
            }

            return I;
        }

        /// Compute the overlap integral between two primitive Gaussians.
        pub fn overlap(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) T {
            const Sij = hermitec(.{self.a[0], other.a[0]}, .{self.alpha, other.alpha}, self.A[0] - other.A[0], 0);
            const Skl = hermitec(.{self.a[1], other.a[1]}, .{self.alpha, other.alpha}, self.A[1] - other.A[1], 0);
            const Smn = hermitec(.{self.a[2], other.a[2]}, .{self.alpha, other.alpha}, self.A[2] - other.A[2], 0);

            return Sij * Skl * Smn * std.math.pow(T, std.math.pi / (self.alpha + other.alpha), 1.5);
        }
    };
}

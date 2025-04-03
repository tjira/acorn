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
            const p = self.alpha + other1.alpha; const q = other2.alpha + other3.alpha; const beta = p * q / (p + q); var j: T = 0;

            const RP: [3]T = .{
                (self.alpha * self.A[0] + other1.alpha * other1.A[0]) / (self.alpha + other1.alpha),
                (self.alpha * self.A[1] + other1.alpha * other1.A[1]) / (self.alpha + other1.alpha),
                (self.alpha * self.A[2] + other1.alpha * other1.A[2]) / (self.alpha + other1.alpha)
            };

            const RQ: [3]T = .{
                (other2.alpha * other2.A[0] + other3.alpha * other3.A[0]) / (other2.alpha + other3.alpha),
                (other2.alpha * other2.A[1] + other3.alpha * other3.A[1]) / (other2.alpha + other3.alpha),
                (other2.alpha * other2.A[2] + other3.alpha * other3.A[2]) / (other2.alpha + other3.alpha)
            };

            const RPQ: [3]T = .{RP[0] - RQ[0], RP[1] - RQ[1], RP[2] - RQ[2]};

            for (0..@as(usize, @intFromFloat(  self.a[0] + other1.a[0])) + 1) |t  | {
                for (0..@as(usize, @intFromFloat(  self.a[1] + other1.a[1])) + 1) |u  | {
                    for (0..@as(usize, @intFromFloat(  self.a[2] + other1.a[2])) + 1) |v  | {
                        for (0..@as(usize, @intFromFloat(other2.a[0] + other3.a[0])) + 1) |tau| {
                            for (0..@as(usize, @intFromFloat(other2.a[1] + other3.a[1])) + 1) |nu | {
                                for (0..@as(usize, @intFromFloat(other2.a[2] + other3.a[2])) + 1) |phi| {

                                    const Eij = hermitec(.{  self.a[0], other1.a[0]}, .{  self.alpha, other1.alpha},   self.A[0] - other1.A[0], @floatFromInt(t  ));
                                    const Ekl = hermitec(.{  self.a[1], other1.a[1]}, .{  self.alpha, other1.alpha},   self.A[1] - other1.A[1], @floatFromInt(u  ));
                                    const Emn = hermitec(.{  self.a[2], other1.a[2]}, .{  self.alpha, other1.alpha},   self.A[2] - other1.A[2], @floatFromInt(v  ));
                                    const Eop = hermitec(.{other2.a[0], other3.a[0]}, .{other2.alpha, other3.alpha}, other2.A[0] - other3.A[0], @floatFromInt(tau));
                                    const Eqr = hermitec(.{other2.a[1], other3.a[1]}, .{other2.alpha, other3.alpha}, other2.A[1] - other3.A[1], @floatFromInt(nu ));
                                    const Est = hermitec(.{other2.a[2], other3.a[2]}, .{other2.alpha, other3.alpha}, other2.A[2] - other3.A[2], @floatFromInt(phi));

                                    const sign = std.math.pow(T, -1, @as(T, @floatFromInt(tau + nu + phi)));

                                    j += sign * Eij * Ekl * Emn * Eop * Eqr * Est * hermitei([3]T{@floatFromInt(t + tau), @floatFromInt(u + nu), @floatFromInt(v + phi)}, RPQ, beta, 0);
                                }
                            }
                        }
                    }
                }
            }

            return 2 * std.math.pow(T, std.math.pi, 2.5) / (p * q * std.math.sqrt(p + q)) * j;
        }

        /// Compute the Hermite Gaussian coefficients in one dimension.
        pub fn hermitec(ij: [2]T, ab: [2]T, Q: T, t: T) T {
            const p = ab[0] + ab[1]; const q = ab[0] * ab[1] / p; var E: T = 0;

            if (ij[0] == 0 and ij[1] == 0 and t == 0) {E += std.math.exp(-q * Q * Q);}

            else if (ij[0] > 0) {

                const E1 = hermitec(.{ij[0] - 1, ij[1]}, ab, Q, t - 1);
                const E2 = hermitec(.{ij[0] - 1, ij[1]}, ab, Q, t    );
                const E3 = hermitec(.{ij[0] - 1, ij[1]}, ab, Q, t + 1);

                E += (1 / (2 * p)) * E1 - (q * Q / ab[0]) * E2 + (t + 1) * E3;
            }

            else if (ij[1] > 0) {

                const E1 = hermitec(.{ij[0], ij[1] - 1}, ab, Q, t - 1);
                const E2 = hermitec(.{ij[0], ij[1] - 1}, ab, Q, t    );
                const E3 = hermitec(.{ij[0], ij[1] - 1}, ab, Q, t + 1);

                E += (1 / (2 * p)) * E1 + (q * Q / ab[1]) * E2 + (t + 1) * E3;
            }

            return E;
        }

        /// Compute the Hermite integrals.
        pub fn hermitei(tuv: [3]T, RPC: [3]T, p: T, n: T) T {
            var I: T = 0;

            if (tuv[0] == 0 and tuv[1] == 0 and tuv[2] == 0) {
                I += std.math.pow(T, -2 * p, n) * mth.boys(p * (RPC[0] * RPC[0] + RPC[1] * RPC[1] + RPC[2] * RPC[2]), n);
            }

            else if (tuv[0] > 0){
                I += (tuv[0] - 1) * hermitei(.{tuv[0] - 2, tuv[1], tuv[2]}, RPC, p, n + 1);
                I +=       RPC[0] * hermitei(.{tuv[0] - 1, tuv[1], tuv[2]}, RPC, p, n + 1);
            }

            else if (tuv[1] > 0) {
                I += (tuv[1] - 1) * hermitei(.{tuv[0], tuv[1] - 2, tuv[2]}, RPC, p, n + 1);
                I +=       RPC[1] * hermitei(.{tuv[0], tuv[1] - 1, tuv[2]}, RPC, p, n + 1);
            }

            else if (tuv[2] > 0) {
                I += (tuv[2] - 1) * hermitei(.{tuv[0], tuv[1], tuv[2] - 2}, RPC, p, n + 1);
                I +=       RPC[2] * hermitei(.{tuv[0], tuv[1], tuv[2] - 1}, RPC, p, n + 1);
            }

            return I;
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
            var n: T = 0;

            for (0..system.atoms.rows) |i| {

                const RPC: [3]T = .{
                    (self.alpha * self.A[0] + other.alpha * other.A[0]) / (self.alpha + other.alpha) - system.coords.at(i, 0),
                    (self.alpha * self.A[1] + other.alpha * other.A[1]) / (self.alpha + other.alpha) - system.coords.at(i, 1),
                    (self.alpha * self.A[2] + other.alpha * other.A[2]) / (self.alpha + other.alpha) - system.coords.at(i, 2)
                };

                for (0..@as(usize, @intFromFloat(self.a[0] + other.a[0])) + 1) |t| {
                    for (0..@as(usize, @intFromFloat(self.a[1] + other.a[1])) + 1) |u| {
                        for (0..@as(usize, @intFromFloat(self.a[2] + other.a[2])) + 1) |v| {

                            const Eij = hermitec(.{self.a[0], other.a[0]}, .{self.alpha, other.alpha}, self.A[0] - other.A[0], @floatFromInt(t));
                            const Ekl = hermitec(.{self.a[1], other.a[1]}, .{self.alpha, other.alpha}, self.A[1] - other.A[1], @floatFromInt(u));
                            const Emn = hermitec(.{self.a[2], other.a[2]}, .{self.alpha, other.alpha}, self.A[2] - other.A[2], @floatFromInt(v));

                            n -= system.atoms.at(i) * Eij * Ekl * Emn * hermitei([3]T{@floatFromInt(t), @floatFromInt(u), @floatFromInt(v)}, RPC, self.alpha + other.alpha, 0);
                        }
                    }
                }
            }

            return 2 * std.math.pi / (self.alpha + other.alpha) * n;
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

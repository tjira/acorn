//! Contracted Gaussian type definitions.

const std = @import("std");

const PrimitiveGaussian = @import("primitivegaussian.zig").PrimitiveGaussian;
const System            = @import("system.zig"           ).System           ;

const dfact = @import("helper.zig").dfact;
const sum   = @import("helper.zig").sum  ;

/// Contracted Gaussian type.
pub fn ContractedGaussian(comptime T: type) type {
    return struct {
        A: [3]T, a: [3]T, c: []T, alpha: []T, l: i8,

        /// Initialize a contracted Gaussian.
        pub fn init(A: [3]T, a: [3]T, c: []const T, alpha: []const T, allocator: std.mem.Allocator) !ContractedGaussian(T) {

            const g = ContractedGaussian(T) {.A = A, .a = a, .c = try allocator.alloc(T, c.len), .alpha = try allocator.alloc(T, alpha.len), .l = @as(i8, @intFromFloat(sum(T, &a)))};

            @memcpy(g.c, c); @memcpy(g.alpha, alpha);

            for (g.c, 0..) |*ci, i| {
                ci.* *= std.math.sqrt(std.math.pow(T, 4 * alpha[i], sum(T, &a)) * std.math.pow(T, 2 * alpha[i] / std.math.pi, 1.5) / dfact(2 * a[0] - 1) / dfact(2 * a[1] - 1) / dfact(2 * a[2] - 1));
            }

            return g;
        }

        /// Deinitialize a contracted Gaussian.
        pub fn deinit(self: ContractedGaussian(T), allocator: std.mem.Allocator) void {
            allocator.free(self.c); allocator.free(self.alpha);
        }

        /// Compute the Coulomb integral between four contracted Gaussians.
        pub fn coulomb(self: ContractedGaussian(T), other1: ContractedGaussian(T), other2: ContractedGaussian(T), other3: ContractedGaussian(T)) T {
            var e: T = 0;

            for (self.c, 0..) |ci, i| for (other1.c, 0..) |cj, j| for (other2.c, 0..) |ck, k| for (other3.c, 0..) |cl, l| {

                const pgi = PrimitiveGaussian(T){.A =   self.A, .a =   self.a, .alpha =   self.alpha[i], .l =   self.l};
                const pgj = PrimitiveGaussian(T){.A = other1.A, .a = other1.a, .alpha = other1.alpha[j], .l = other1.l};
                const pgk = PrimitiveGaussian(T){.A = other2.A, .a = other2.a, .alpha = other2.alpha[k], .l = other2.l};
                const pgl = PrimitiveGaussian(T){.A = other3.A, .a = other3.a, .alpha = other3.alpha[l], .l = other3.l};

                e += ci * cj * ck * cl * pgi.coulomb(pgj, pgk, pgl);
            };

            return e;
        }

        /// Compute the kinetic integral between two contracted Gaussians.
        pub fn kinetic(self: ContractedGaussian(T), other: ContractedGaussian(T)) T {
            var t: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A =  self.A, .a =  self.a, .alpha =  self.alpha[i], .l =  self.l};
                const pgj = PrimitiveGaussian(T){.A = other.A, .a = other.a, .alpha = other.alpha[j], .l = other.l};

                t += ci * cj * pgi.kinetic(pgj);
            };

            return t;
        }

        /// Compute the nuclear integral between two contracted Gaussians.
        pub fn nuclear(self: ContractedGaussian(T), other: ContractedGaussian(T), system: System(T)) T {
            var v: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A =  self.A, .a =  self.a, .alpha =  self.alpha[i], .l =  self.l};
                const pgj = PrimitiveGaussian(T){.A = other.A, .a = other.a, .alpha = other.alpha[j], .l = other.l};

                v += ci * cj * pgi.nuclear(pgj, system);
            };

            return v;
        }

        /// Compute the overlap integral between two contracted Gaussians.
        pub fn overlap(self: ContractedGaussian(T), other: ContractedGaussian(T)) T {
            var s: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A =  self.A, .a =  self.a, .alpha =  self.alpha[i], .l =  self.l};
                const pgj = PrimitiveGaussian(T){.A = other.A, .a = other.a, .alpha = other.alpha[j], .l = other.l};

                s += ci * cj * pgi.overlap(pgj);
            };

            return s;
        }
    };
}

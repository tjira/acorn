//! Module for checking and generating prime numbers

const std = @import("std");

const inp = @import("input.zig" );

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The main function for prime checking and generation
pub fn run(comptime T: type, opt: inp.PrimeOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    _ = allocator; var timer = try std.time.Timer.start();

    if (std.mem.eql(u8, opt.mode, "basic")) {

        if (opt.generate == null and print) {

            const divisor = smallest_divisor(T, opt.number);

            try std.io.getStdOut().writer().print("\nNUMBER: {d}, SMALLEST DIVISOR: {d}, PRIME: {any}\n", .{opt.number, divisor, divisor == 1 and opt.number > 1});
        }

        if (opt.generate != null and print) {

            try std.io.getStdOut().writer().print("\nGENERATED PRIMES\n", .{});

            var n = nextprime(T, opt.number); var i: T = 0;

            while (i < opt.generate.?.count) : ({i += 1; n = nextprime(T, n);}) {
                if ((i + 1) % opt.generate.?.log_interval == 0 or i == 0) try std.io.getStdOut().writer().print("{d}: {d} ({s})\n", .{i + 1, n, std.fmt.fmtDuration(timer.read())});
            }
        }
    }

    else if (std.mem.eql(u8, opt.mode, "mersenne")) {

        if (opt.generate == null and print) {
            try std.io.getStdOut().writer().print("\nEXPONENT: {d}, MERSENNE PRIME: {any}\n", .{opt.number, ismersenne(T, opt.number)});
        }

        if (opt.generate != null and print) {

            try std.io.getStdOut().writer().print("\nGENERATED MERSENNE PRIME EXPONENTS\n", .{});

            var p = nextmersenne(T, opt.number); var i: T = 0;

            while (i < opt.generate.?.count) : ({i += 1; p = nextmersenne(T, p);}) {
                if ((i + 1) % opt.generate.?.log_interval == 0 or i == 0) try std.io.getStdOut().writer().print("{d}: {d} ({s})\n", .{i + 1, p, std.fmt.fmtDuration(timer.read())});
            }
        }
    }

    else return error.InvalidPrimeMode;

}

/// The Lucas-Lehmer test for mersenne primes
pub fn ismersenne(comptime T: type, p: T) bool {
    if (smallest_divisor(T, p) != 1) return false;

    const M = std.math.pow(T, 2, p) - 1; var s: T = 4;

    var i: T = 0; while (i < p - 2) : (i += 1) s = (s * s - 2) % M;

    return s == 0 or p == 2;
}

/// The find the next mersenne prime
pub fn nextmersenne(comptime T: type, p: T) T {
    if (p < 2) return 2;

    var np = p + 2; if (p % 2 == 0) np -= 1;

    while (!ismersenne(T, np)) np += 2;

    return np;
}

/// The find the next prime
pub fn nextprime(comptime T: type, n: T) T {
    if (n < 2) return 2;

    var nn = n + 2; if (n % 2 == 0) nn -= 1;

    while (smallest_divisor(T, nn) != 1) nn += 2;

    return nn;
}

/// The find the smallest divisor of a number
pub fn smallest_divisor(comptime T: type, n: T) T {
    if (n % 2 == 0) return 2;

    var d: T = 3; while (d * d <= n) : (d += 2) if (n % d == 0) return d;

    return 1;
}

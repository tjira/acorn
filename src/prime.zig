//! Module for checking and generating prime numbers

const std = @import("std");

const inp = @import("input.zig" );

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The main function for prime checking and generation
pub fn run(comptime T: type, opt: inp.PrimeOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    var timer = try std.time.Timer.start(); var p = opt.number - 1;

    if (std.mem.eql(u8, opt.mode, "basic")) {

        if (opt.generate == null and print) {

            try std.io.getStdOut().writer().print("\nNUMBER: {d}, PRIME: {any}\n", .{opt.number, isprime(p)});
        }

        if (opt.generate != null and print) {

            try std.io.getStdOut().writer().print("\nPRIME NUMBERS\n", .{});

            for (0..opt.generate.?.count) |i| {

                p = nextprime(p);

                if ((i + 1) % opt.generate.?.log_interval == 0 or i == 0) {
                    try std.io.getStdOut().writer().print("{d}: {d} ({s})\n", .{i + 1, p, std.fmt.fmtDuration(timer.read())});
                }
            }
        }
    }

    else if (std.mem.eql(u8, opt.mode, "mersenne")) {

        if (opt.generate == null and print) {
            try std.io.getStdOut().writer().print("\nEXPONENT: {d}, MERSENNE PRIME: {any}\n", .{opt.number, try ismersenne(allocator, p)});
        }

        if (opt.generate != null and print) {

            try std.io.getStdOut().writer().print("\nMERSENNE PRIME EXPONENTS\n", .{});

            for (0..opt.generate.?.count) |i| {

                p = try nextmersenne(allocator, p);

                if ((i + 1) % opt.generate.?.log_interval == 0 or i == 0) {
                    try std.io.getStdOut().writer().print("{d}: {d} ({s})\n", .{i + 1, p, std.fmt.fmtDuration(timer.read())});
                }
            }
        }
    }

    else return error.InvalidPrimeMode;
}

/// The Lucas-Lehmer test for Mersenne primes
pub fn ismersenne(allocator: std.mem.Allocator, p: anytype) !bool {
    if (!isprime(p)) return false;

    var M = try std.math.big.int.Managed.initSet(allocator, 2); defer M.deinit();
    var s = try std.math.big.int.Managed.initSet(allocator, 4); defer s.deinit();
    var q = try std.math.big.int.Managed.init   (allocator   ); defer q.deinit();

    try std.math.big.int.Managed.shiftLeft(&M, &M, @as(u32, @intCast(p - 1)));
    try std.math.big.int.Managed.addScalar(&M, &M, @as(i32, @intCast(0 - 1)));

    for (0..@as(usize, @intCast(p - 2))) |_| {

        try std.math.big.int.Managed.mul      (&s, &s, &s);
        try std.math.big.int.Managed.addScalar(&s, &s, -2);

        try std.math.big.int.Managed.divFloor(&q, &s, &s, &M);
    }

    return std.math.big.int.Managed.eqlZero(s) or p == 2;
}

/// The find the smallest divisor of a number
pub fn isprime(p: anytype) bool {
    if (p < 2 or (p > 2 and p % 2 == 0)) return false;

    var d: @TypeOf(p) = 3;

    while (d * d <= p) : (d += 2) if (p % d == 0) return false;

    return true;
}

/// The find the next mersenne prime
pub fn nextmersenne(allocator: std.mem.Allocator, p: anytype) !@TypeOf(p) {
    if (p < 2) return 2; var np = p + 1 + p % 2;

    while (!try ismersenne(allocator, np)) np += 2;

    return np;
}

/// The find the next prime
pub fn nextprime(p: anytype) @TypeOf(p) {
    if (p < 2) return 2; var np = p + 1 + p % 2;

    while (!isprime(np)) np += 2;

    return np;
}

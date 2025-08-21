//! Module for checking and generating prime numbers.

const std = @import("std");

const hlp = @import("helper.zig");
const inp = @import("input.zig" );

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The main function for prime checking and generation.
pub fn run(comptime T: type, opt: inp.PrimeOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    var timer = try std.time.Timer.start(); var p = opt.number - 1; var P: Matrix(T) = undefined;

    if (opt.generate != null and opt.generate.?.output != null) {
        P = try Matrix(T).init(opt.generate.?.count, 2, allocator);
    }

    if (std.mem.eql(u8, opt.mode, "basic")) {

        if (opt.generate == null and print) {
            try hlp.print(std.fs.File.stdout(), "\nNUMBER: {d}, PRIME: {any}\n", .{opt.number, isPrime(p)});
        }

        if (opt.generate != null and print) {

            try hlp.print(std.fs.File.stdout(), "\nPRIME NUMBERS\n", .{});

            for (0..opt.generate.?.count) |i| {

                p = nextPrime(p);

                if (opt.generate.?.output != null) {
                    P.ptr(i, 0).* = @as(T, @intCast(i)) + 1; P.ptr(i, 1).* = p;
                }

                if ((i + 1) % opt.generate.?.log_interval == 0 or i == 0) {
                    try hlp.print(std.fs.File.stdout(), "{d}: {d} ({D})\n", .{i + 1, p, timer.read()});
                }
            }
        }
    }

    else if (std.mem.eql(u8, opt.mode, "mersenne")) {

        if (opt.generate == null and print) {
            try hlp.print(std.fs.File.stdout(), "\nEXPONENT: {d}, MERSENNE PRIME: {any}\n", .{opt.number, try isMersenne(allocator, p)});
        }

        if (opt.generate != null and print) {

            try hlp.print(std.fs.File.stdout(), "\nMERSENNE PRIME EXPONENTS\n", .{});

            for (0..opt.generate.?.count) |i| {

                p = try nextMersenne(allocator, p);

                if (opt.generate.?.output != null) {
                    P.ptr(i, 0).* = @as(T, @intCast(i)) + 1; P.ptr(i, 1).* = p;
                }

                if ((i + 1) % opt.generate.?.log_interval == 0 or i == 0) {
                    try hlp.print(std.fs.File.stdout(), "{d}: {d} ({D})\n", .{i + 1, p, timer.read()});
                }
            }
        }
    }

    else return error.InvalidPrimeMode;

    if (opt.generate != null and opt.generate.?.output != null) {
        try P.write(opt.generate.?.output.?); P.deinit();
    }
}

/// The Lucas-Lehmer test for Mersenne primes.
pub fn isMersenne(allocator: std.mem.Allocator, p: anytype) !bool {
    if (!isPrime(p)) return false;

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

/// Check if the provided number is prime.
pub fn isPrime(p: anytype) bool {
    if (p < 2 or (p > 2 and p % 2 == 0)) return false;

    var d: @TypeOf(p) = 3;

    while (d * d <= p) : (d += 2) if (p % d == 0) return false;

    return true;
}

/// The find the next Mersenne prime.
pub fn nextMersenne(allocator: std.mem.Allocator, p: anytype) !@TypeOf(p) {
    if (p < 2) return 2; var np = p + 1 + p % 2;

    while (!try isMersenne(allocator, np)) np += 2;

    return np;
}

/// The find the next prime.
pub fn nextPrime(p: anytype) @TypeOf(p) {
    if (p < 2) return 2; var np = p + 1 + p % 2;

    while (!isPrime(np)) np += 2;

    return np;
}

//! File for math functions.

const std = @import("std");

const gsl = @import("gsl.zig");

const asfloat = @import("helper.zig").asfloat;

/// Boys function.
pub fn boys(comptime T: type, a: T, x: T) T {
    return gsl.hyp1f1(a + 0.5, a + 1.5, -x) / (2 * a + 1);
}

/// Calculate the combination number.
pub fn comb(n: anytype, k: @TypeOf(n)) @TypeOf(n) {
    var nck: @TypeOf(n) = 1;

    for (k + 1..n + 1) |i| nck *= i;
    for (2..n - k + 1) |i| nck /= i;

    return nck;
}

/// Generate all combinations of n elements from an array.
pub fn combinations(comptime T: type, array: []const T, n: usize, allocator: std.mem.Allocator) !std.ArrayList([]T) {
    var result = std.ArrayList([]T).init(allocator); var current = std.ArrayList(usize).init(allocator); defer current.deinit();

    const Backtrack = struct { fn get(res: *std.ArrayList([]T), arr: []const T, m: usize, curr: *std.ArrayList(usize), start: usize, alloc: std.mem.Allocator) !void {
        if (curr.items.len == m) {

            var last = try alloc.alloc(T, m);

            for (curr.items, 0..) |index, i| {
                last[i] = arr[index];
            }

            try res.append(last); return;
        }

        for (start..arr.len) |i| {
            try curr.append(i); try get(res, arr, m, curr, i + 1, alloc); _ = curr.pop();
        }
    }};

    try Backtrack.get(&result, array, n, &current, 0, allocator); return result;
}

/// Return the maximum of two numbers
pub fn max(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (a > b) a else b;
}

/// Return the minimum of two numbers
pub fn min(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (a < b) a else b;
}

/// Calculate the product of an array.
pub fn prod(comptime T: type, v: []const T) T {
    var result: T = 1;

    for (v) |value| result *= value;

    return result;
}

/// Sign of a number.
pub fn sgn(a: anytype) @TypeOf(a) {
    return if (a < 0) -1 else 1;
}

/// Calculate the sum of an array.
pub fn sum(comptime T: type, v: []const T) T {
    var result: T = 0;

    for (v) |value| result += value;

    return result;
}

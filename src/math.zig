//! File for math functions.

const std = @import("std");

const asfloat = @import("helper.zig").asfloat;

/// Calculate the combination number.
pub fn comb(n: anytype, k: @TypeOf(n)) @TypeOf(n, k) {
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

/// Gamma function using the Lanczos approximation.
pub fn gamma(comptime T: type, x: std.math.Complex(T)) std.math.Complex(T) {
    const g: T = 7; const p = [9]T{
           0.999999999999809930e+0,
          676.52036812188510000e+0,
        -1259.13921672240280000e+0,
          771.32342877765313000e+0,
         -176.61502916214059000e+0,
           12.50734327868690500e+0,
           -0.13857109526572012e+0,
            9.98436957801957160e-6,
            1.50563273514931160e-7
    };

    if (x.re < 0.5) {
        return std.math.Complex(T).init(std.math.pi, 0).div(gamma(T, std.math.Complex(T).init(1, 0).sub(x)).mul(std.math.complex.sin(x.mul(std.math.Complex(T).init(std.math.pi, 0)))));
    }

    var a = std.math.Complex(T).init(p[0], 0);

    for (1..g + 2) |i| a = a.add(std.math.Complex(T).init(p[i], 0).div(x.add(std.math.Complex(T).init(asfloat(T, i) - 1, 0))));

    const t = x.add(std.math.Complex(T).init(g - 0.5, 0));

    const b = std.math.Complex(T).init(std.math.sqrt(2 * std.math.pi), 0);

    return a.mul(b).mul(std.math.complex.pow(t, x.sub(std.math.Complex(T).init(0.5, 0)))).mul(std.math.complex.exp(t.neg()));
}

/// Return the Heaviside function
pub fn h(x: anytype) @TypeOf(x) {
    return if (x >= 0) 1.0 else 0;
}

/// Return the maximum of two numbers
pub fn max(a: anytype, b: @TypeOf(a)) @TypeOf(a, b) {
    return if (a > b) a else b;
}

/// Calculates the mean of an array.
pub fn mean(comptime T: type, v: []const T) !T {
    if (v.len == 0) return error.InvalidInputForMean;

    var result: T = 0;

    for (v) |value| {
        result += value;
    }

    return result / asfloat(T, v.len);
}

/// Return the minimum of two numbers
pub fn min(a: anytype, b: @TypeOf(a)) @TypeOf(a, b) {
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

/// Calculates the standard deviation of an array.
pub fn stddev(comptime T: type, v: []const T) !T {
    if (v.len < 2) return error.InvalidInputForStddev;

    const mean_value = try mean(T, v); var sum_squared_diff: T = 0;

    for (v) |value| {
        sum_squared_diff += value * value + mean_value * mean_value - 2 * value * mean_value;
    }

    return std.math.sqrt(sum_squared_diff / asfloat(T, v.len - 1));
}

/// Calculate the sum of an array.
pub fn sum(comptime T: type, v: []const T) T {
    var result: T = 0;

    for (v) |value| result += value;

    return result;
}

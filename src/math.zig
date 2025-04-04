//! File for math functions.

const std = @import("std");

const asfloat = @import("helper.zig").asfloat;

/// Boys function with zero n.
pub fn boys(comptime T: type, x: T, a: T) T {
    return if (x > 0) 0.5 * gamma(T, a + 0.5) * gammainc(T, x, a + 0.5) / std.math.pow(T, x, a + 0.5) else 1 / (2 * a + 1);
}

/// Calculate the combination number.
pub fn comb(n: anytype, k: @TypeOf(n)) @TypeOf(n) {
    var nck: @TypeOf(n) = 1;

    for (k + 1..n + 1) |i| nck *= i;
    for (2..n - k + 1) |i| nck /= i;

    return nck;
}

/// Calculate the double factorial of a number.
pub fn dfact(n: anytype) @TypeOf(n) {
    if (n == -1 or n == 0 or n == 1) return 1;

    if (n == 2) return 2;

    return n * dfact(n - 2);
}

/// Calculate the factorial of a number.
pub fn fact(n: usize) @TypeOf(n) {
    var result: @TypeOf(n) = 1;

    for (2..n + 1) |i| result *= i;

    return result;
}

/// Gamma function using the Lanczos approximation.
pub fn gamma(comptime T: type, x: T) T {
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

    if (x < 0.5) return std.math.pi / (std.math.sin(std.math.pi * x) * gamma(T, 1 - x));
    
    var a: T = p[0];

    for (1..g + 2) |i| a += p[i] / (x + asfloat(T, i) - 1);

    const t: T = x + g - 0.5;

    return std.math.sqrt(2 * std.math.pi) * std.math.pow(T, t, x - 0.5) * std.math.exp(-t) * a;
}

/// Regularized lower incomplete gamma function.
pub fn gammainc(comptime T: type, x: T, a: T) T {
    const eps: T = 1e-14;

    if (x < a + 1) {

        var b = 1 / a; var c = 1 / a; var i: T = 1;

        while (true) : (i += 1) {
            c *= x / (a + i); b += c; if (c < eps) break;
        }

        return b * std.math.exp(-x + a * std.math.log(T, std.math.e, x)) / gamma(T, a);
    }

    else {

        var b = x + 1 - a; var c: T = 1e10; var d = 1 / b; var e = d; var i: T = 1;

        while (true) : (i += 1) {
            const f = i * (a - i); b += 2; d = f * d + b; c = b + f / c; d = 1 / d; const g = d * c; e *= g; if (@abs(g - 1) < eps) break;
        }

        return 1 - std.math.exp(-x + a * std.math.log(T, std.math.e, x) - std.math.log(T, std.math.e, gamma(T, a))) * e;
    }
}

/// Return the maximum of two numbers
pub fn max(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (a > b) a else b;
}

/// Return the mean of a vector.
pub fn mean(comptime T: type, v: []const T) T {
    return sum(T, v) / asfloat(T, v.len);
}

/// Return the minimum of two numbers
pub fn min(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (a < b) a else b;
}

/// Calculate the power of a number.
pub fn powi(a: anytype, b: u32) @TypeOf(a) {
    var result: @TypeOf(a) = 1;

    for (0..b) |_| result *= a;

    return result;
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

/// Return the standard deviation of a vector.
pub fn stdev(comptime T: type, v: []const T) T {
    var result: T = 0;

    for (v) |value| result += (value - mean(T, v)) * (value - mean(T, v));

    return std.math.sqrt(result / asfloat(T, v.len));
}

/// Calculate the sum of an array.
pub fn sum(comptime T: type, v: []const T) T {
    var result: T = 0;

    for (v) |value| result += value;

    return result;
}

//! Module with some helper functions.

const std = @import("std");

/// Casts an integer to a float.
pub fn asfloat(comptime T: type, a: anytype) T {
    return @as(T, @floatFromInt(a));
}

/// Reverse the bits of a number.
pub fn bitrev(value: anytype, count: u6) @TypeOf(value) {
    var result: @TypeOf(value) = 0; var i: u6 = 0; while (i < count) : (i += 1) {result |= ((value >> @intCast(i)) & 1) << @intCast(count - 1 - i);} return result;
}

/// Calculate the combination number.
pub fn c(n: anytype, k: @TypeOf(n)) @TypeOf(n) {
    var nck: @TypeOf(n) = 1; for (k + 1..n + 1) |i| {nck *= i;} for (2..n - k + 1) |i| {nck /= i;} return nck;
}

/// Return the maximum of two numbers
pub fn max(comptime T: type, a: T, b: T) T {
    return if (a > b) a else b;
}

/// Return the minimum of two numbers
pub fn min(comptime T: type, a: T, b: T) T {
    return if (a < b) a else b;
}

/// Calculate the product of an array.
pub fn prod(comptime T: type, v: []const T) T {
    var result: T = 1; for (v) |value| result *= value; return result;
}

/// Remove the carriage return from a string. Fucking windows.
pub fn uncr(string: []const u8) []const u8 {
    return if (string[string.len - 1] == 13) string[0..string.len - 1] else string;
}

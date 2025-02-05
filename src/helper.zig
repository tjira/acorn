//! Module with some helper functions.

const std = @import("std");

/// Absolute value of a number.
pub fn abs(a: anytype) @TypeOf(a) {
    return if (a < 0) -a else a;
}

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

/// Check if a vector contains a value.
pub fn contains(comptime T: type, v: []const T, value: T) bool {
    for (v) |element| {if (element == value) return true;} return false;
}

// Calculate the factorial of a number.
pub fn fact(n: usize) @TypeOf(n) {
    var result: @TypeOf(n) = 1; for (2..n + 1) |i| {result *= i;} return result;
}

/// Check if the passed type is a struct.
pub fn istruct(comptime T: type) bool {
    return @typeInfo(T) == .Struct;
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

/// Sign of a number.
pub fn sgn(a: anytype) @TypeOf(a) {
    return if (a < 0) -1 else 1;
}

/// Calculate the sum of an array.
pub fn sum(comptime T: type, v: []const T) T {
    var result: T = 0; for (v) |value| result += value; return result;
}

/// Remove the carriage return from a string. Fucking windows.
pub fn uncr(string: []const u8) []const u8 {
    return if (string[string.len - 1] == 13) string[0..string.len - 1] else string;
}

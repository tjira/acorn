//! Module with some helper functions.

const std = @import("std");

/// Casts an integer to a float.
pub fn asfloat(comptime T: type, a: anytype) T {
    return @as(T, @floatFromInt(a));
}

/// Check if a vector contains a value.
pub fn contains(comptime T: type, v: []const T, value: T) bool {
    for (v) |element| {if (element == value) return true;} return false;
}

/// Check if the passed type is a struct.
pub fn istruct(comptime T: type) bool {
    return @typeInfo(T) == .@"struct";
}

/// Remove the carriage return from a string. Fucking windows.
pub fn uncr(string: []const u8) []const u8 {
    return if (string[string.len - 1] == 13) string[0..string.len - 1] else string;
}

//! Module with some helper functions.

const std = @import("std");

/// Casts an integer to a float.
pub fn asfloat(comptime T: type, a: anytype) T {
    return @as(T, @floatFromInt(a));
}

/// Reverse the bits of a number.
pub fn bitrev(value: anytype, count: u6) @TypeOf(value) {
    var result: @TypeOf(value) = 0; var i: u6 = 0;

    while (i < count) : (i += 1) {
        result |= ((value >> @intCast(i)) & 1) << @intCast(count - 1 - i);
    }

    return result;
}

/// Check if a vector contains a value.
pub fn contains(comptime T: type, v: []const T, value: T) bool {
    for (v) |element| {if (element == value) return true;} return false;
}

/// Check if the passed type is a struct.
pub fn istruct(comptime T: type) bool {
    return @typeInfo(T) == .@"struct";
}

/// Return the maximum of two numbers.
pub fn max(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (a > b) a else b;
}

/// Return the formatted memory size of a type and count.
pub fn memFormat(comptime T: type, count: usize, allocator: std.mem.Allocator) ![]const u8 {
    const buffer = try allocator.alloc(u8, 9); const bytes = @sizeOf(T) * count;

    const digits = std.math.log10(bytes) + 1;

    const units = switch (digits) {
        1...3 => " B", 4...6 => "kB", 7...9 => "MB", 10...12 => "GB", 13...15 => "TB", else => unreachable
    };

    const size = asfloat(f64, bytes) / switch (digits) {
        1...3 => @as(f64, 1e0), 4...6 => 1e3, 7...9 => 1e6, 10...12 => 1e9, 13...15 => 1e12, else => unreachable
    };

    return try std.fmt.bufPrint(buffer, "{d:6.2} {s}", .{size, units});
}

/// Remove the carriage return from a string. Fucking windows.
pub fn uncr(string: []const u8) []const u8 {
    return if (string[string.len - 1] == 13) string[0..string.len - 1] else string;
}

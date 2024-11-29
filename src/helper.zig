const std = @import("std");

pub fn asfloat(comptime T: type, a: anytype) T {
    return @as(T, @floatFromInt(a));
}

pub fn bitrev(value: anytype, count: u6) @TypeOf(value) {
    var result: @TypeOf(value) = 0; var i: u6 = 0; while (i < count) : (i += 1) {result |= ((value >> i) & 1) << (count - 1 - i);} return result;
}

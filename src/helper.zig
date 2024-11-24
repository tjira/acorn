const std = @import("std");

pub fn asfloat(comptime T: type, a: anytype) T {
    return @as(T, @floatFromInt(a));
}

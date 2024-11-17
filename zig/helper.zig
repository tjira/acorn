const std = @import("std");

pub fn asfloat(comptime T: type, a: anytype) T {
    return @as(T, @floatFromInt(a));
}

pub fn swap(a: anytype, b: anytype) void {
    const c = a.*; a.* = b.*; b.* = c;
}

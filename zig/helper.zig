const std = @import("std");

pub fn swap(a: anytype, b: anytype) void {
    const c = a.*; a.* = b.*; b.* = c;
}

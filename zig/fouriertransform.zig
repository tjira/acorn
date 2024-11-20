const std = @import("std");

const asfloat = @import("helper.zig").asfloat;

pub fn dft(comptime T: type, out: []std.math.Complex(T), in: []const std.math.Complex(T), shape: []const usize, factor: T) void {
    const size = shape[0]; var exp = std.math.Complex(T).init(0, 0); for (0..size) |i| out[i] = std.math.Complex(T).init(0, 0);

    for (0..size) |i| for (0..size) |j| {
        exp.im = factor * 2 * std.math.pi * asfloat(T, i * j) / asfloat(T, size); out[i] = out[i].add(std.math.complex.exp(exp).mul(in[j]));
    };

    if (factor < 0) for (0..size) |i| {out[i] = out[i].div(std.math.Complex(T).init(asfloat(T, size), 0));};
}

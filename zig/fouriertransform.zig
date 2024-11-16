const std = @import("std");

pub fn dft(comptime T: type, output: []std.math.Complex(T), input: []const std.math.Complex(T), shape: []const usize) void {
    const length = shape[0];

    for (0..length) |i| output[i] = std.math.Complex(T).init(0, 0);

    for (0..length) |i| {
        for (0..length) |j| {

            const exponent = std.math.Complex(T).init(0, -2 * std.math.pi * @as(T, @floatFromInt(i * j)) / @as(T, @floatFromInt(length)));

            output[i] = output[i].add(std.math.complex.exp(exponent).mul(input[j]));
        }
    }
}

const std = @import("std");

const matrix = @import("matrix.zig");

pub fn main() !void {
    const samples = 10;

    const matrix_mm   = try matrix.matrix_mm  (samples); std.debug.print("MATRIX MULTIPLICATION: {d:.3} +- {d:.3}\n", .{matrix_mm[0],   matrix_mm[1]  });
    const matrix_eigh = try matrix.matrix_eigh(samples); std.debug.print("MATRIX EIGENPROBLEM:   {d:.3} +- {d:.3}\n", .{matrix_eigh[0], matrix_eigh[1]});
}

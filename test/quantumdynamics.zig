const std = @import("std");

const mat = @import("acorn").mat;
const qdn = @import("acorn").qdn;

const Matrix = @import("acorn").Matrix;

const allocator = std.testing.allocator;

test "qdyn_default" {
    const opt = qdn.QuantumDynamicsOptions(f64){};

    const output = try qdn.run(f64, opt, false, allocator); defer output.deinit();

    try std.testing.expect(@abs(output.Ekin + output.Epot - 1.12490708654035) < 1e-8);
}

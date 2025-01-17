const std = @import("std");

const cdn = @import("acorn").cdn;
const vec = @import("acorn").vec;

const Vector = @import("acorn").Vector;

const allocator = std.testing.allocator;

test "cdyn_default" {
    const opt = cdn.ClassicalDynamicsOptions(f64){};

    const output = try cdn.run(f64, opt, false, allocator); defer output.deinit();

    try std.testing.expect(@abs(output.Ekin + output.Epot - 1.09553019324015) < 1e-8);
}

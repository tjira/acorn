const std = @import("std");

const matrix             = @import("matrix.zig"           );
const quantum_dynamics   = @import("quantumdynamics.zig"  );
const classical_dynamics = @import("classicaldynamics.zig");
const tensor             = @import("tensor.zig"           );
const vector             = @import("vector.zig"           );

comptime {
    _ = classical_dynamics; _ = matrix; _ = quantum_dynamics; _ = tensor; _ = vector;
}

test {
    std.testing.refAllDecls(@This());
}

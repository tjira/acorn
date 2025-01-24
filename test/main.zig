const std = @import("std");

const matrix = @import("matrix.zig"           );
const tensor = @import("tensor.zig"           );
const vector = @import("vector.zig"           );

comptime {
    _ = matrix;
    _ = tensor;
    _ = vector;
}

test {
    std.testing.refAllDecls(@This());
}

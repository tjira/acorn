const std = @import("std");

const mat = @import("matrix.zig"           );
const ten = @import("tensor.zig"           );
const vec = @import("vector.zig"           );

comptime {
    _ = mat;
    _ = ten;
    _ = vec;
}

test {
    std.testing.refAllDecls(@This());
}

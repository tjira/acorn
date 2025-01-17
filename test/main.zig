const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const mat = @import("matrix.zig"           );
const qdn = @import("quantumdynamics.zig"  );
const ten = @import("tensor.zig"           );
const vec = @import("vector.zig"           );

comptime {
    _ = cdn;
    _ = mat;
    _ = qdn;
    _ = ten;
    _ = vec;
}

test {
    std.testing.refAllDecls(@This());
}

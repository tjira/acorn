const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const ftr = @import("fouriertransform.zig" );
const mat = @import("matrix.zig"           );
const qdn = @import("quantumdynamics.zig"  );
const vec = @import("vector.zig"           );

comptime {
    _ = cdn; _ = ftr; _ = mat; _ = qdn; _ = vec;
}

test {
    std.testing.refAllDecls(@This());
}

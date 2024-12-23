const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const qdn = @import("quantumdynamics.zig"  );

comptime {
    _ = cdn; _ = qdn;
}

test {
    std.testing.refAllDecls(@This());
}

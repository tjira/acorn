const std = @import("std");

const math               = @import("math.zig"             );
const matrix             = @import("matrix.zig"           );
const quantum_dynamics   = @import("quantumdynamics.zig"  );
const classical_dynamics = @import("classicaldynamics.zig");
const tensor             = @import("tensor.zig"           );
const vector             = @import("vector.zig"           );

pub fn expect(val: anytype, ref: @TypeOf(val)) !void {
    std.testing.expect(@abs(val - ref) < 1e-12) catch |err| {
        std.debug.print("VALUE: {d:.14}, REFERENCE: {d:.14}\n", .{val, ref}); return err;
    };
}

pub fn log(val: anytype, ref: @TypeOf(val)) void {
    std.debug.print("VALUE: {d:.14}, REFERENCE: {d:.14}\n", .{val, ref});
}

comptime {
    _ = classical_dynamics; _ = math; _ = matrix; _ = quantum_dynamics; _ = tensor; _ = vector;
}

test {
    std.testing.refAllDecls(@This());
}

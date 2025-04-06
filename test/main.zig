const std = @import("std");

const classical_dynamics        = @import("classicaldynamics.zig"       );
const hartree_fock              = @import("hartreefock.zig"             );
const math                      = @import("math.zig"                    );
const matrix                    = @import("matrix.zig"                  );
const moller_plesset            = @import("mollerplesset.zig"           );
const configuration_interaction = @import("configurationinteraction.zig");
const quantum_dynamics          = @import("quantumdynamics.zig"         );
const tensor                    = @import("tensor.zig"                  );
const vector                    = @import("vector.zig"                  );

pub fn expect(val: anytype, ref: @TypeOf(val)) !void {
    std.testing.expect(@abs(val - ref) < 1e-12) catch |err| {
        std.debug.print("VALUE: {d:.14}, REFERENCE: {d:.14}\n", .{val, ref}); return err;
    };
}

pub fn log(val: anytype, ref: @TypeOf(val)) void {
    std.debug.print("VALUE: {d:.14}, REFERENCE: {d:.14}\n", .{val, ref});
}

comptime {
    _ = classical_dynamics; _ = configuration_interaction; _ = hartree_fock; _ = math; _ = matrix; _ = moller_plesset; _ = quantum_dynamics; _ = tensor; _ = vector;
}

test {
    std.testing.refAllDecls(@This());
}

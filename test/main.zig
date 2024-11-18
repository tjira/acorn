const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const mat = @import("matrix.zig"           );
const mpt = @import("modelpotential.zig"   );
const qdn = @import("quantumdynamics.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.testing.allocator; const tol = 1e-14;

test "matrixFill" {
    try std.testing.expect(true);
}

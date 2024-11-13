const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn ModelPotential(comptime T: type) type {
    return struct {
        pub fn tully(r: Vector(T)) !Matrix(T) {
            const U = try Matrix(T).zero(1, 1, r.data.allocator);
            try U.set(0, 0, 0.5 * try r.at(0) * try r.at(0));
            return U;
        }
    };
}

//! File to calculate properties from a quantum mechanical calculation

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Vector = @import("vector.zig").Vector;

/// Returns the harmonic frequencies of the system.
pub fn freq(comptime T: type, system: System(T), H: Matrix(T), allocator: std.mem.Allocator) !Vector(T) {
    const f = try Vector(T).init(H.rows, allocator); f.fill(0);

    _ = system;

    return f;
}

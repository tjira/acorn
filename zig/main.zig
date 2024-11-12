const std = @import("std");

const ClassicalDynamics = @import("classicaldynamics.zig").ClassicalDynamics;
const Matrix            = @import("matrix.zig"           ).Matrix           ;
const Vector            = @import("vector.zig"           ).Vector           ;

const allocator = std.heap.page_allocator;

pub fn potential_tully(r: Vector(f64)) !Matrix(f64) {
    const U = try Matrix(f64).zero(1, 1, allocator); // define the zero potential
    try U.set(0, 0, 0.5 * try r.at(0) * try r.at(0)); // set the value at (0, 0) of the potential
    return U; // return the potential
}

pub fn main() !void {
    const r0m = try Matrix(f64).init(1, 1, &[_]f64{2}, allocator);
    const v0m = try Matrix(f64).init(1, 1, &[_]f64{0}, allocator);
    const a0m = try Matrix(f64).init(1, 1, &[_]f64{0}, allocator);
    const s0v = try Vector(u8).init(1, &[_]u8{0}, allocator);
    const classicalDynamics = ClassicalDynamics{.iterations = 10000, .mass = 1, .timeStep = 0.01};
    try classicalDynamics.run(potential_tully, r0m, v0m, a0m, s0v);
}

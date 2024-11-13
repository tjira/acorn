const std = @import("std");

const ClassicalDynamics = @import("classicaldynamics.zig").ClassicalDynamics;
const Matrix            = @import("matrix.zig"           ).Matrix           ;
const ModelPotential    = @import("modelpotential.zig"   ).ModelPotential   ;
const Vector            = @import("vector.zig"           ).Vector           ;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    const r0 = try Matrix(f64).init(1, 1, &[_]f64{2}, allocator);
    const v0 = try Matrix(f64).init(1, 1, &[_]f64{0}, allocator);
    const a0 = try Matrix(f64).init(1, 1, &[_]f64{0}, allocator);
    const s0 = try Vector(u8).init(1, &[_]u8{0}, allocator);
    const classical_dynamics = ClassicalDynamics(f64){.iterations = 10000, .mass = 1, .time_step = 0.01};
    try classical_dynamics.run(ModelPotential(f64).tully, r0, v0, a0, s0);

    const A = try Matrix(f64).init(3, 2, &[_]f64{1, 2, 3, 4, 5, 6}, allocator);
    try A.write("A.mat");
}

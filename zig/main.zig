const std = @import("std");

const ClassicalDynamics = @import("classicaldynamics.zig").ClassicalDynamics;
const Matrix            = @import("matrix.zig"           ).Matrix           ;
const ModelPotential    = @import("modelpotential.zig"   ).ModelPotential   ;
const Vector            = @import("vector.zig"           ).Vector           ;

const evaluateModelPotential = @import("modelpotential.zig").evaluateModelPotential;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    const potential = ModelPotential(f64).tully1;
    const r0 = try Matrix(f64).randNorm(1, 1, -10, 0.5, 1, allocator);
    const p0 = try Matrix(f64).randNorm(1, 1,  12, 1.0, 1, allocator);
    const a0 = try Matrix(f64).zero    (1, 1,              allocator);
    const s0 = try Vector(u8 ).constant(1, 1,              allocator);

    const classical_dynamics = ClassicalDynamics(f64){.iterations = 10, .mass = 2000, .time_step = 0.01, .adiabatic = false, .seed = 1};
    try classical_dynamics.run(potential, r0, try p0.divScalar(classical_dynamics.mass), a0, s0);

    try (try evaluateModelPotential(f64, potential, -16, 16, 1024, allocator)).write("POTENTIAL.mat");
}

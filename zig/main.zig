const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const mpt = @import("modelpotential.zig"   );
const mat = @import("matrix.zig"           );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const global_allocator = std.heap.page_allocator;

pub fn main() !void {
    const cdyn_opt = cdn.ClassicalDynamicsOptions(f64){
        .adiabatic = false,
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 1,
        .seed = 1,
        .ic = .{
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .state = 1,
            .mass = 2000
        },
        .potential = mpt.doubleState1D_1,
    };

    try cdn.run(f64, cdyn_opt, global_allocator);
}

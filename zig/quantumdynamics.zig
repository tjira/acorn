const std = @import("std");

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn QuantumDynamicsOptions(comptime T: type) type {
    return struct {
        const Grid = struct {
            limits: []const T, points: u32
        };
        const InitialConditions = struct {
            position: []const T, momentum: []const T, state: u32, mass: T
        };
        const LogIntervals = struct {
            iteration: u32
        };

        iterations: u32,
        time_step: f64,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals, potential: mpt.PotentialType(T),
    };
}

pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    _ = opt;
    _ = allocator;
}

const std = @import("std");

const mpt = @import("modelpotential.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        const InitialConditions = struct {
            position_mean: []const T, position_std: []const T, momentum_mean: []const T, momentum_std: []const T, state: u8, mass: T
        };

        adiabatic: bool,
        iterations: u32,
        time_step: T,
        trajectories: u32,
        seed: u64,

        ic: InitialConditions, potential: mpt.PotentialType(T),
    };
}
pub fn run(comptime T: type, opt: ClassicalDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var prng = std.Random.DefaultPrng.init(opt.seed); const rand = prng.random();

    for (0..opt.trajectories) |_| {

        var r = try Vector(f64).init(1, allocator); defer r.deinit();
        var p = try Vector(f64).init(1, allocator); defer p.deinit();
        var a = try Vector(f64).init(1, allocator); defer a.deinit();
        var s = opt.ic.state; a.fill(0);

        s = s;

        for (0..r.rows) |i| r.ptr(i).* = opt.ic.position_mean[i] + opt.ic.position_std[i] * rand.floatNorm(T);
        for (0..p.rows) |i| p.ptr(i).* = opt.ic.momentum_mean[i] + opt.ic.momentum_std[i] * rand.floatNorm(T);

        var v = try p.clone(); defer v.deinit(); v.ptr(0).* /= opt.ic.mass;

        var U = try Matrix(f64).init(2, 2, allocator); defer U.deinit();
        var F = try Vector(f64).init(r.rows, allocator); defer F.deinit();

        for (0..opt.iterations) |_| {

            try mpt.potentialAndForce(T, &U, &F, opt.potential, r, s);

            for (0..r.rows) |k| {
                const ap = try a.clone();
                a.ptr(k).* = F.at(k) / opt.ic.mass;
                v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step;
            }

            std.debug.print("{any}\n", .{r.data});
        }
    }
}

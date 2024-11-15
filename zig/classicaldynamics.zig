const std = @import("std");

const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        const InitialConditions = struct {
            position_mean: []const T, position_std: []const T, momentum_mean: []const T, momentum_std: []const T, state: u32, mass: T
        };
        const LogIntervals = struct {
            trajectory: u32, iteration: u32
        };

        adiabatic: bool,
        derivative_step: T,
        iterations: u32,
        seed: u32,
        time_step: T,
        trajectories: u32,

        ic: InitialConditions, li: LogIntervals, potential: mpt.PotentialType(T),
    };
}

pub fn run(comptime T: type, opt: ClassicalDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var prng = std.Random.DefaultPrng.init(opt.seed); const rand = prng.random(); const nstate = try mpt.states(T, opt.potential);

    var pop = try Matrix(T).init(opt.iterations, nstate + 1, allocator); defer pop.deinit(); pop.fill(0);

    for (0..opt.iterations) |i| pop.ptr(i, 0).* = @as(T, @floatFromInt(i + 1)) * opt.time_step;

    for (0..opt.trajectories) |i| {

        var r  = try Vector(T).init(1, allocator); defer  r.deinit();
        var p  = try Vector(T).init(1, allocator); defer  p.deinit();
        var a  = try Vector(T).init(1, allocator); defer  a.deinit();
        var ap = try Vector(T).init(1, allocator); defer ap.deinit();
        var s = opt.ic.state; var sp = s; a.fill(0);

        for (0..r.rows) |j| r.ptr(j).* = opt.ic.position_mean[j] + opt.ic.position_std[j] * rand.floatNorm(T);
        for (0..p.rows) |j| p.ptr(j).* = opt.ic.momentum_mean[j] + opt.ic.momentum_std[j] * rand.floatNorm(T);

        var v = try p.clone(); defer v.deinit(); v.ptr(0).* /= opt.ic.mass; var Ekin: T = 0; var Epot: T = 0;

        var U = try Matrix(T).init(nstate, nstate, allocator); defer U.deinit();
        var F = try Vector(T).init(r.rows,         allocator); defer F.deinit();

        var U3 = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer U3[0].deinit(); defer U3[1].deinit(); defer U3[2].deinit();

        for (0..opt.iterations) |j| {

            if (j == 0 or sp != s) try mpt.evaluate(T, &U, &F, opt.potential, r, s, opt.derivative_step);

            @memcpy(ap.data, a.data); sp = s; pop.ptr(j, 1 + s).* += 1;

            for (0..r.rows) |k| {
                a.ptr(k).* = F.at(k) / opt.ic.mass;
                v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step;
            }

            try mpt.evaluate(T, &U, &F, opt.potential, r, s, opt.derivative_step); @memcpy(U3[j % 3].data, U.data);

            Ekin = 0; for (v.data) |e| {Ekin += e * e;} Ekin *= 0.5 * opt.ic.mass; Epot = U.at(s, s);

            if (j > 2) s = try landauZener(T, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, rand);

            if ((i == 0 or (i + 1) % opt.li.trajectory == 0) and (j == 0 or (j + 1) % opt.li.iteration == 0)) {
                std.debug.print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:4} {d:12.6}\n", .{i + 1, j + 1, Ekin, Epot, Ekin + Epot, sp, r.at(0)});
            }
        }
    }

    for (0..opt.iterations) |i| {for (0..try mpt.states(T, opt.potential)) |j| pop.ptr(i, j + 1).* /= opt.trajectories;}

    for (0..try mpt.states(T, opt.potential)) |i| {
        std.debug.print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i + 1)});
    }

    try mat.write(T, "POPULATION.mat", pop);
}

fn landauZener(comptime T: type, U3: []const Matrix(T), s: u32, time_step: T, rand: std.Random) !u32 {
    var dU = try std.ArrayList(T).initCapacity(U3[0].allocator, U3[0].rows); defer dU.deinit();

    for (0..U3[0].rows) |i| try dU.append((U3[0].at(i, i) - U3[1].at(i, i)) / time_step);

    const p = 1 - @exp(-2 * std.math.pi * std.math.pow(T, U3[0].at(0, 1), 2) / @abs(dU.items[1] - dU.items[0]));

    if ((U3[0].at(1, 1) - U3[0].at(0, 0)) * (U3[1].at(1, 1) - U3[1].at(0, 0)) < 0 and rand.float(T) < p) return (s + 1) % 2;

    return s;
}

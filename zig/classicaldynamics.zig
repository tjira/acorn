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

    var pop = try Matrix(T).init(opt.iterations, nstate, allocator); defer pop.deinit(); pop.fill(0);

    var r  = try Vector(T).init(1, allocator); defer  r.deinit();
    var p  = try Vector(T).init(1, allocator); defer  p.deinit();
    var v  = try Vector(T).init(1, allocator); defer  v.deinit();
    var a  = try Vector(T).init(1, allocator); defer  a.deinit();
    var rt = try Vector(T).init(1, allocator); defer rt.deinit();
    var ap = try Vector(T).init(1, allocator); defer ap.deinit();

    var U = try Matrix(T).init(nstate, nstate, allocator); defer U.deinit();
    var A = try Matrix(T).init(nstate, nstate, allocator); defer A.deinit();
    var C = try Matrix(T).init(nstate, nstate, allocator); defer C.deinit();

    var U3 = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer U3[0].deinit(); defer U3[1].deinit(); defer U3[2].deinit();

    for (0..opt.trajectories) |i| {

        for (0..r.rows) |j| r.ptr(j).* = opt.ic.position_mean[j] + opt.ic.position_std[j] * rand.floatNorm(T);
        for (0..p.rows) |j| p.ptr(j).* = opt.ic.momentum_mean[j] + opt.ic.momentum_std[j] * rand.floatNorm(T);

        v = try p.clone(); v.ptr(0).* /= opt.ic.mass; a.fill(0); var s = opt.ic.state; var sp = s; var Ekin: T = 0; var Epot: T = 0;

        for (0..opt.iterations) |j| {

            @memcpy(ap.data, a.data); sp = s; pop.ptr(j, s).* += 1;

            for (0..r.rows) |k| {

                rt.ptr(k).* = r.at(k) + opt.derivative_step; opt.potential(T, &U, rt); const Up = U.at(s, s);
                rt.ptr(k).* = r.at(k) - opt.derivative_step; opt.potential(T, &U, rt); const Um = U.at(s, s);

                a.ptr(k).* = -0.5 * (Up - Um) / opt.derivative_step / opt.ic.mass;
                v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step;
            }

            opt.potential(T, &U, r); @memcpy(U3[j % 3].data, U.data); Ekin = 0; for (v.data) |e| {Ekin += e * e;} Ekin *= 0.5 * opt.ic.mass; Epot = U.at(s, s);

            if (j > 2) s = try landauZener(T, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, opt.adiabatic, rand);

            if ((i == 0 or (i + 1) % opt.li.trajectory == 0) and (j == 0 or (j + 1) % opt.li.iteration == 0)) {
                std.debug.print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:4} {d:12.6}\n", .{i + 1, j + 1, Ekin, Epot, Ekin + Epot, sp, r.at(0)});
            }
        }
    }

    for (0..opt.iterations) |i| {for (0..try mpt.states(T, opt.potential)) |j| pop.ptr(i, j).* /= opt.trajectories;}

    for (0..try mpt.states(T, opt.potential)) |i| {
        std.debug.print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});
    }

    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * opt.iterations);

    var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); time.hjoin(&pop_t, pop); try pop_t.write("POPULATION.mat"); pop_t.deinit();
}

fn landauZener(comptime T: type, U3: []const Matrix(T), s: u32, time_step: T, adiabatic: bool, rand: std.Random) !u32 {
    var sn = s; var pm: T = 0; const rn = rand.float(T); _=adiabatic;

    for (0..U3[0].rows) |i| if (i != s and (U3[0].at(i, i) - U3[0].at(s, s)) * (U3[1].at(i, i) - U3[1].at(s, s)) < 0) {

        const di = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const ds = (U3[0].at(s, s) - U3[1].at(s, s)) / time_step;

        const p = 1 - std.math.exp(-2 * std.math.pi * std.math.pow(T, U3[0].at(i, s), 2) / @abs(di - ds));

        if (rn < p and p > pm) {sn = @intCast(i); pm = p;}
    };

    return sn;
}

const std = @import("std"); const gsl = @cImport(@cInclude("gsl/gsl_eigen.h"));

const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        const InitialConditions = struct {
            position_mean: []const T, position_std: []const T, momentum_mean: []const T, momentum_std: []const T, state: u32, mass: T
        };
        const LogIntervals = struct {
            trajectory: u32, iteration: u32
        };
        const Write = struct {
            population: ?[]const u8
        };

        adiabatic: bool,
        derivative_step: T,
        iterations: u32,
        seed: u32,
        time_step: T,
        trajectories: u32,

        initial_conditions: InitialConditions, log_intervals: LogIntervals, write: Write, potential: []const u8
    };
}

pub fn run(comptime T: type, opt: ClassicalDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var prng = std.Random.DefaultPrng.init(opt.seed); const rand = prng.random();

    var pop = try Matrix(T).init(opt.iterations, std.math.pow(u32, mpt.states(opt.potential), 1), allocator); defer pop.deinit(); pop.fill(0);
    var tdc = try Matrix(T).init(opt.iterations, std.math.pow(u32, mpt.states(opt.potential), 2), allocator); defer tdc.deinit(); tdc.fill(0);

    {
        const GSLW = gsl.gsl_eigen_symmv_alloc(mpt.states(opt.potential)); defer gsl.gsl_eigen_symmv_free(GSLW);

        var r  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  r.deinit();
        var p  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  p.deinit();
        var v  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  v.deinit();
        var a  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  a.deinit();
        var rt = try Vector(T).init(mpt.dims(opt.potential), allocator); defer rt.deinit();
        var ap = try Vector(T).init(mpt.dims(opt.potential), allocator); defer ap.deinit();

        var U    = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer    U.deinit();
        var UA   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer   UA.deinit();
        var UC   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer   UC.deinit();
        var UCS  = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  UCS.deinit();
        var UCSA = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer UCSA.deinit();
        var UCSC = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer UCSC.deinit();
        var TDC  = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  TDC.deinit();

        var U3  = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer  U3[0].deinit(); defer  U3[1].deinit(); defer U3[2].deinit();
        var UC2 = [2]Matrix(T){try U.clone(), try U.clone()               }; defer UC2[0].deinit(); defer UC2[1].deinit()                      ;

        for (0..opt.trajectories) |i| {

            for (0..r.rows) |j| r.ptr(j).* = opt.initial_conditions.position_mean[j] + opt.initial_conditions.position_std[j] * rand.floatNorm(T);
            for (0..p.rows) |j| p.ptr(j).* = opt.initial_conditions.momentum_mean[j] + opt.initial_conditions.momentum_std[j] * rand.floatNorm(T);

            @memcpy(v.data, p.data); v.ptr(0).* /= opt.initial_conditions.mass; a.fill(0); var s = opt.initial_conditions.state; var sp = s; var Ekin: T = 0; var Epot: T = 0;

            for (0..opt.iterations) |j| {

                @memcpy(ap.data, a.data); sp = s; pop.ptr(j, s).* += 1;

                for (0..r.rows) |k| {

                    rt.ptr(k).* = r.at(k) + opt.derivative_step; mpt.eval(T, &U, opt.potential, rt); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, GSLW); @memcpy(U.data, UA.data);} const Up = U.at(s, s);
                    rt.ptr(k).* = r.at(k) - opt.derivative_step; mpt.eval(T, &U, opt.potential, rt); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, GSLW); @memcpy(U.data, UA.data);} const Um = U.at(s, s);

                    a.ptr(k).* = -0.5 * (Up - Um) / opt.derivative_step / opt.initial_conditions.mass;
                    v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                    r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step;
                }

                mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {

                    mat.eigh(T, &UA, &UC, U, GSLW); @memcpy(U.data, UA.data); @memcpy(UC2[j % 2].data, UC.data);

                    if (j > 1) for (0..UC.cols) |k| {
                        var overlap: T = 0; for (0..UC.rows) |l| {overlap += UC2[j % 2].at(l, k) * UC2[(j - 1) % 2].at(l, k);} if (overlap < 0) for (0..UC.rows) |l| {UC2[j % 2].ptr(l, k).* *= -1;};
                    };

                    TDC.fill(0); if (j > 1) derivativeCoupling(T, &TDC, &[_]Matrix(T){UC2[j % 2], UC2[(j - 1) % 2]}, opt.time_step, &UCS, &UCSA, &UCSC, GSLW);

                    for (0..TDC.data.len) |k| tdc.ptr(j, k).* = UCS.data[k];
                }

                Ekin = 0; for (v.data) |e| {Ekin += e * e;} Ekin *= 0.5 * opt.initial_conditions.mass; Epot = U.at(s, s);

                if ((i == 0 or (i + 1) % opt.log_intervals.trajectory == 0) and (j == 0 or (j + 1) % opt.log_intervals.iteration == 0)) {
                    try std.io.getStdOut().writer().print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:4} {d:12.6}\n", .{i + 1, j + 1, Ekin, Epot, Ekin + Epot, s, r.at(0)});
                }

                @memcpy(U3[j % 3].data, U.data);

                if (j > 1) s = try landauZener(T, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, opt.adiabatic, rand);
                // if (j > 0) s = try fewestSwitches(T, TDC, s, opt.time_step, rand);

                if (s != sp and Ekin < U.at(s, s) - U.at(sp, sp)) s = sp;

                if (sp != s) {for (0..v.rows) |k| v.ptr(k).* *= std.math.sqrt((Ekin - U.at(s, s) + U.at(sp, sp)) / Ekin);}
            }
        }
    }

    for (0..opt.iterations) |i| {for (0..mpt.states(opt.potential)) |j| pop.ptr(i, j).* /= asfloat(T, opt.trajectories);}

    for (0..mpt.states(opt.potential)) |i| {
        try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});
    }

    try writeResults(T, opt, pop, tdc, allocator);
}

fn derivativeCoupling(comptime T: type, TDC: *Matrix(T), UC2: []const Matrix(T), time_step: T, UCS: *Matrix(T), UCSA: *Matrix(T), UCSC: *Matrix(T), GSLW: *gsl.gsl_eigen_symmv_workspace) void {
    UCS.fill(0);

    for (0..UCS.rows) |i| for (0..UCS.cols) |j| for (0..UCS.rows) |k| {
        UCS.ptr(i, j).* += UC2[1].at(k, i) * UC2[0].at(k, j);
    };

    mat.eigh(T, UCSA, UCSC, UCS.*, GSLW);

    for (0..UCSA.rows) |i| UCSA.ptr(i, i).* = std.math.log10(UCSA.at(i, i)) / std.math.log10e;

    mat.mm(T, TDC, UCSC.*, UCSA.*); UCSC.transpose(); mat.mm(T, UCSA, TDC.*, UCSC.*); @memcpy(TDC.data, UCSA.data);

    for (0..TDC.rows) |i| for (0..TDC.cols) |j| {TDC.ptr(i, j).* /= time_step;};
}

fn fewestSwitches(comptime T: type, TDC: Matrix(T), s: u32, time_step: T, rand: std.Random) !u32 {
    _ = TDC; _ = time_step; _ = rand;

    return s;
}

fn landauZener(comptime T: type, U3: []const Matrix(T), s: u32, time_step: T, adiabatic: bool, rand: std.Random) !u32 {
    var sn = s; var pm: T = 0; var rn: T = undefined; var sampled = false;

    if (!adiabatic) for (0..U3[0].rows) |i| if (i != s) {

        if ((U3[0].at(i, i) - U3[0].at(s, s)) * (U3[1].at(i, i) - U3[1].at(s, s)) > 0) continue;

        const di = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const ds = (U3[0].at(s, s) - U3[1].at(s, s)) / time_step;

        const p = 1 - std.math.exp(-2 * std.math.pi * std.math.pow(T, U3[0].at(i, s), 2) / @abs(di - ds));

        if (!sampled) {rn = rand.float(T); sampled = true;} if (rn < p and p > pm) {sn = @intCast(i); pm = p;}
    };

    if (adiabatic) for (0..U3[0].rows) |i| if (i != s) {

        const di0 = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const ds0 = (U3[0].at(s, s) - U3[1].at(s, s)) / time_step;
        const di1 = (U3[1].at(i, i) - U3[2].at(i, i)) / time_step; const ds1 = (U3[1].at(s, s) - U3[2].at(s, s)) / time_step;

        if ((di0 - ds0) > 0 or (di1 - ds1) < 0) continue;

        const ddi = (di0 - di1) / time_step; const dds = (ds0 - ds1) / time_step;

        const p = std.math.exp(-0.5 * std.math.pi * std.math.sqrt(std.math.pow(T, U3[0].at(i, i) - U3[0].at(s, s), 3) / (ddi - dds)));

        if (!sampled) {rn = rand.float(T); sampled = true;} if (rn < p and p > pm) {sn = @intCast(i); pm = p;}
    };

    return sn;
}

fn writeResults(comptime T: type, opt: ClassicalDynamicsOptions(T), pop: Matrix(T), tdc: Matrix(T), allocator: std.mem.Allocator) !void {
    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * asfloat(T, opt.iterations));

    if (opt.write.population) |path| {
        var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); mat.hjoin(T, &pop_t, time, pop); try pop_t.write(path); pop_t.deinit();
    }

    var tdc_t = try Matrix(T).init(opt.iterations, tdc.cols + 1, allocator); mat.hjoin(T, &tdc_t, time, tdc); try tdc_t.write("TDC.mat"); tdc_t.deinit();
}

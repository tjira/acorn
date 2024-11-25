const std = @import("std"); const Complex = std.math.Complex; const gsl_eigen = @cImport(@cInclude("gsl/gsl_eigen.h"));

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
            kinetic_energy_mean: ?[]const u8,
            momentum_mean: ?[]const u8,
            population_mean: ?[]const u8,
            position_mean: ?[]const u8,
            potential_energy_mean: ?[]const u8,
            total_energy_mean: ?[]const u8
        };

        adiabatic: bool,
        derivative_step: T,
        iterations: u32,
        potential: []const u8,
        seed: u32,
        time_step: T,
        trajectories: u32,
        type: []const u8,

        initial_conditions: InitialConditions, log_intervals: LogIntervals, write: Write, 
    };
}

pub fn run(comptime T: type, opt: ClassicalDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var prng = std.Random.DefaultPrng.init(opt.seed); const rand = prng.random();

    const fssh = std.mem.eql(u8, opt.type, "fssh"); const kfssh = std.mem.eql(u8, opt.type, "kfssh"); const lzsh = std.mem.eql(u8, opt.type, "lzsh");

    var pop      = try Matrix(T).init(opt.iterations, mpt.states(opt.potential), allocator); defer      pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations, 1                        , allocator); defer     ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations, 1                        , allocator); defer     epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations, 1                        , allocator); defer     etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer momentum.deinit(); momentum.fill(0);

    {
        const GSLW = gsl_eigen.gsl_eigen_symmv_alloc(mpt.states(opt.potential)); defer gsl_eigen.gsl_eigen_symmv_free(GSLW);

        var r  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  r.deinit();
        var p  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  p.deinit();
        var v  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  v.deinit();
        var a  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  a.deinit();
        var ap = try Vector(T).init(mpt.dims(opt.potential), allocator); defer ap.deinit();

        var U    = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer   U.deinit();
        var UA   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  UA.deinit();
        var UC   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  UC.deinit();
        var UT   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  UT.deinit();
        var UCS  = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer UCS.deinit();
        var TDC  = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer TDC.deinit();

        var C = try Vector(Complex(T)).init(mpt.states(opt.potential), allocator); defer C.deinit(); 

        var K1 = try Vector(Complex(T)).init(C.rows, allocator); defer K1.deinit();
        var K2 = try Vector(Complex(T)).init(C.rows, allocator); defer K2.deinit();
        var K3 = try Vector(Complex(T)).init(C.rows, allocator); defer K3.deinit();
        var K4 = try Vector(Complex(T)).init(C.rows, allocator); defer K4.deinit();

        var U3  = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer  U3[0].deinit(); defer  U3[1].deinit(); defer U3[2].deinit();
        var UC2 = [2]Matrix(T){try U.clone(), try U.clone()               }; defer UC2[0].deinit(); defer UC2[1].deinit()                      ;

        try std.io.getStdOut().writer().print("\n{s:6} {s:6} {s:12} {s:12} {s:12} {s:5}", .{"TRAJ", "ITER", "EKIN", "EPOT", "ETOT", "STATE"});

        if (r.rows > 1) for (0..r.rows - 1) |_| {try std.io.getStdOut().writer().print(" " ** 14, .{});}; try std.io.getStdOut().writer().print(" {s:14}",   .{"POSITION"  });
        if (v.rows > 1) for (0..v.rows - 1) |_| {try std.io.getStdOut().writer().print(" " ** 14, .{});}; try std.io.getStdOut().writer().print(" {s:14}\n", .{"MOMENTUM"  });

        for (0..opt.trajectories) |i| {

            for (0..r.rows) |j| r.ptr(j).* = opt.initial_conditions.position_mean[j] + opt.initial_conditions.position_std[j] * rand.floatNorm(T);
            for (0..p.rows) |j| p.ptr(j).* = opt.initial_conditions.momentum_mean[j] + opt.initial_conditions.momentum_std[j] * rand.floatNorm(T);

            C.fill(Complex(T).init(0, 0)); C.ptr(opt.initial_conditions.state).* = Complex(T).init(1, 0);

            @memcpy(v.data, p.data); v.ptr(0).* /= opt.initial_conditions.mass; a.fill(0); var s = opt.initial_conditions.state; var sp = s; var Ekin: T = 0; var Epot: T = 0;

            for (0..opt.iterations) |j| {

                @memcpy(ap.data, a.data); sp = s; if (opt.write.population_mean != null) pop.ptr(j, s).* += 1;

                for (0..r.rows) |k| {

                    r.ptr(k).* += 1 * opt.derivative_step; mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, &UT, GSLW); @memcpy(U.data, UA.data);} const Up = U.at(s, s);
                    r.ptr(k).* -= 2 * opt.derivative_step; mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, &UT, GSLW); @memcpy(U.data, UA.data);} const Um = U.at(s, s);

                    a.ptr(k).* = -0.5 * (Up - Um) / opt.derivative_step / opt.initial_conditions.mass;
                    v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                    r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step + opt.derivative_step;
                }

                mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {

                    mat.eigh(T, &UA, &UC, U, &UT, GSLW); @memcpy(U.data, UA.data); @memcpy(UC2[j % 2].data, UC.data);

                    if (j > 0) for (0..UC.cols) |k| {
                        var overlap: T = 0; for (0..UC.rows) |l| {overlap += UC2[j % 2].at(l, k) * UC2[(j - 1) % 2].at(l, k);} if (overlap < 0) for (0..UC.rows) |l| {UC2[j % 2].ptr(l, k).* *= -1;};
                    };
                }

                @memcpy(U3[j % 3].data, U.data); Ekin = 0; for (v.data) |e| {Ekin += e * e;} Ekin *= 0.5 * opt.initial_conditions.mass; Epot = U.at(s, s);

                if (opt.adiabatic and  fssh and j > 0) derivativeCouplingNumeric(T, &TDC, &UCS, &[_]Matrix(T){UC2[j % 2], UC2[(j - 1) % 2]},                opt.time_step);
                if (opt.adiabatic and kfssh and j > 1) derivativeCouplingBaeckan(T, &TDC,       &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, opt.time_step);

                if (opt.write.kinetic_energy_mean   != null) ekin.ptr(j, 0).* += Ekin                                                            ;
                if (opt.write.potential_energy_mean != null) epot.ptr(j, 0).* += Epot                                                            ;
                if (opt.write.total_energy_mean     != null) etot.ptr(j, 0).* += Ekin + Epot                                                     ;
                if (opt.write.position_mean         != null) for (0..r.rows) |k| {position.ptr(j, k).* += r.at(k);}                              ;
                if (opt.write.momentum_mean         != null) for (0..v.rows) |k| {momentum.ptr(j, k).* += v.at(k) * opt.initial_conditions.mass;};

                if ((i == 0 or (i + 1) % opt.log_intervals.trajectory == 0) and (j == 0 or (j + 1) % opt.log_intervals.iteration == 0)) {
                    try printIteration(T, @intCast(i), @intCast(j), Ekin, Epot, Ekin + Epot, s, r, v, opt.initial_conditions.mass);
                }

                if ((lzsh         ) and j > 1) s = try landauZener(T, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, opt.adiabatic, rand);
                if ((fssh or kfssh) and j > 1) s = try fewestSwitches(T, &C, U, TDC, s, opt.time_step, rand, &K1, &K2, &K3, &K4);

                if (s != sp and Ekin < U.at(s, s) - U.at(sp, sp)) s = sp;

                if (sp != s) {for (0..v.rows) |k| v.ptr(k).* *= std.math.sqrt((Ekin - U.at(s, s) + U.at(sp, sp)) / Ekin);}
            }
        }
    }

    for (0..opt.iterations) |i| {for (0..pop.cols     ) |j|     pop.ptr(i, j).*  /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..ekin.cols    ) |j|     ekin.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..epot.cols    ) |j|     epot.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..etot.cols    ) |j|     etot.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..position.cols) |j| position.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..momentum.cols) |j| momentum.ptr(i, j).* /= asfloat(T, opt.trajectories);}

    for (0..mpt.states(opt.potential)) |i| {
        try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});
    }

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, allocator);
}

fn derivativeCouplingBaeckan(comptime T: type, TDC: *Matrix(T), U3: []const Matrix(T), time_step: T) void {
    for (0..TDC.rows) |i| for (i + 1..TDC.cols) |j| {

        const di0 = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const dj0 = (U3[0].at(j, j) - U3[1].at(j, j)) / time_step;
        const di1 = (U3[1].at(i, i) - U3[2].at(i, i)) / time_step; const dj1 = (U3[1].at(j, j) - U3[2].at(j, j)) / time_step;

        const ddi = (di0 - di1) / time_step; const ddj = (dj0 - dj1) / time_step; const arg = (ddi - ddj) / (U3[0].at(i, i) - U3[0].at(j, j));

        if (arg > 0) {TDC.ptr(i, j).* = 0.5 * std.math.sqrt(arg);} TDC.ptr(j, i).* = -TDC.at(i, j);
    };
}

fn derivativeCouplingNumeric(comptime T: type, TDC: *Matrix(T), UCS: *Matrix(T), UC2: []const Matrix(T), time_step: T) void {
    UCS.fill(0);

    for (0..UCS.rows) |i| for (0..UCS.cols) |j| for (0..UCS.rows) |k| {
        UCS.ptr(i, j).* += UC2[1].at(k, i) * UC2[0].at(k, j);
    };

    for (0..TDC.rows) |i| for (0..TDC.cols) |j| {
        TDC.ptr(i, j).* = (UCS.at(i, j) - UCS.at(j, i)) / (2 * time_step);
    };
}

fn fewestSwitches(comptime T: type, C: *Vector(Complex(T)), U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, rand: std.Random, K1: @TypeOf(C), K2: @TypeOf(C), K3: @TypeOf(C), K4: @TypeOf(C)) !u32 {
    const iters = 10; var ns = s;

    const Function = struct { fn get (K: *Vector(Complex(T)), FC: Vector(Complex(T)), FU: Matrix(T), FTDC: Matrix(T)) void {
        for (0..FC.rows) |i| {K.ptr(i).* = FC.at(i).mul(Complex(T).init(FU.at(i, i), 0)).mulbyi().neg();} for (0..FC.rows) |i| for (0..FC.rows) |j| {
            K.ptr(i).* = K.at(i).sub(FC.at(j).mul(Complex(T).init(FTDC.at(i, j), 0)));
        };
    }};

    for (0..iters) |_| {

        K1.fill(Complex(T).init(0, 0)); K2.fill(Complex(T).init(0, 0)); K3.fill(Complex(T).init(0, 0)); K4.fill(Complex(T).init(0, 0));

        Function.get(K1, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).add(K1.at(j).mul(Complex(T).init(time_step / 2 / iters, 0)));
        }

        Function.get(K2, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).sub(K1.at(j).mul(Complex(T).init(time_step / 2 / iters, 0)));
            C.ptr(j).* = C.at(j).add(K2.at(j).mul(Complex(T).init(time_step / 2 / iters, 0)));
        }

        Function.get(K3, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).sub(K2.at(j).mul(Complex(T).init(time_step / 2 / iters, 0)));
            C.ptr(j).* = C.at(j).add(K3.at(j).mul(Complex(T).init(time_step / 1 / iters, 0)));
        }

        Function.get(K4, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).sub(K3.at(j).mul(Complex(T).init(time_step / iters, 0)));
        }

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).add(K1.at(j).add(K2.at(j).mul(Complex(T).init(2, 0))).add(K3.at(j).mul(Complex(T).init(2, 0))).add(K4.at(j)).mul(Complex(T).init(time_step / 6 / iters, 0)));
        }

        const rn = rand.float(T); var pm: T = 0; for (0..C.rows) |j| if (j != ns) {
            const p = 2 * TDC.at(ns, j) * C.at(j).mul(C.at(ns).conjugate()).re / std.math.pow(T, C.at(ns).magnitude(), 2) * time_step / iters; if (rn < p and p > pm) {ns = @intCast(j); pm = p;}
        };
    }

    return ns;
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

fn printIteration(comptime T: type, i: u32, j: u32, Ekin: T, Epot: T, Etot: T, s: u32, r: Vector(T), v: Vector(T), mass: T) !void {
    try std.io.getStdOut().writer().print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:5} [", .{i + 1, j + 1, Ekin, Epot, Etot, s});

    for (0..r.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:12.4}", .{if (k == 0) "" else ", ", r.at(k)});
    }

    try std.io.getStdOut().writer().print("] [", .{});

    for (0..v.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:12.4}", .{if (k == 0) "" else ", ", v.at(k) * mass});
    }

    try std.io.getStdOut().writer().print("]\n", .{});
}

fn writeResults(comptime T: type, opt: ClassicalDynamicsOptions(T), pop: Matrix(T), ekin: Matrix(T), epot: Matrix(T), etot: Matrix(T), position: Matrix(T), momentum: Matrix(T), allocator: std.mem.Allocator) !void {
    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * asfloat(T, opt.iterations));

    if (opt.write.kinetic_energy_mean) |path| {
        var ekin_t = try Matrix(T).init(opt.iterations, ekin.cols + 1, allocator); mat.hjoin(T, &ekin_t, time, ekin); try ekin_t.write(path); ekin_t.deinit();
    }

    if (opt.write.population_mean) |path| {
        var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); mat.hjoin(T, &pop_t, time, pop); try pop_t.write(path); pop_t.deinit();
    }

    if (opt.write.potential_energy_mean) |path| {
        var epot_t = try Matrix(T).init(opt.iterations, epot.cols + 1, allocator); mat.hjoin(T, &epot_t, time, epot); try epot_t.write(path); epot_t.deinit();
    }

    if (opt.write.total_energy_mean) |path| {
        var etot_t = try Matrix(T).init(opt.iterations, etot.cols + 1, allocator); mat.hjoin(T, &etot_t, time, etot); try etot_t.write(path); etot_t.deinit();
    }

    if (opt.write.position_mean) |path| {
        var position_t = try Matrix(T).init(opt.iterations, position.cols + 1, allocator); mat.hjoin(T, &position_t, time, position); try position_t.write(path); position_t.deinit();
    }

    if (opt.write.momentum_mean) |path| {
        var momentum_t = try Matrix(T).init(opt.iterations, momentum.cols + 1, allocator); mat.hjoin(T, &momentum_t, time, momentum); try momentum_t.write(path); momentum_t.deinit();
    }
}

const std = @import("std"); const Complex = std.math.Complex;

const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        const InitialConditions = struct {
            position_mean: []const T = &[_]f64{-10}, position_std: []const T = &[_]f64{0.5}, momentum_mean: []const T = &[_]f64{15}, momentum_std: []const T = &[_]f64{1}, state: u32 = 1, mass: T = 2000
        };
        const LogIntervals = struct {
            trajectory: u32 = 0, iteration: u32 = 0
        };
        const Write = struct {
            fssh_coefficient_mean: ?[]const u8 = null,
            kinetic_energy_mean: ?[]const u8 = null,
            momentum_mean: ?[]const u8 = null,
            population_mean: ?[]const u8 = null,
            position_mean: ?[]const u8 = null,
            potential_energy_mean: ?[]const u8 = null,
            total_energy_mean: ?[]const u8 = null
        };

        adiabatic: bool = true,
        derivative_step: T = 0.001,
        iterations: u32 = 3000,
        seed: u32 = 1,
        time_step: T = 1,
        trajectories: u32 = 100,

        initial_conditions: InitialConditions = .{}, log_intervals: LogIntervals = .{}, write: Write = .{}, potential: []const u8 = "tully1D_1", type: []const u8 = "fssh"
    };
}

pub fn ClassicalDynamicsOutput(comptime T: type) type {
    return struct {
        pop: Vector(T),

        pub fn init(nstate: usize, allocator: std.mem.Allocator) !ClassicalDynamicsOutput(T) {
            return ClassicalDynamicsOutput(T){
                .pop = try Vector(T).init(nstate, allocator)
            };
        }
        pub fn deinit(self: ClassicalDynamicsOutput(T)) void {
            self.pop.deinit();
        }
    };
}

pub fn run(comptime T: type, opt: ClassicalDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !ClassicalDynamicsOutput(T) {
    var output = try ClassicalDynamicsOutput(T).init(mpt.states(opt.potential), allocator);

    var prng_jump = std.Random.DefaultPrng.init(opt.seed); const rand_jump = prng_jump.random();
    var prng_traj = std.Random.DefaultPrng.init(opt.seed); const rand_traj = prng_traj.random();
    var prng_bloc = std.Random.DefaultPrng.init(opt.seed); const rand_bloc = prng_bloc.random();

    const fssh = std.mem.eql(u8, opt.type, "fssh"); const kfssh = std.mem.eql(u8, opt.type, "kfssh");
    const mash = std.mem.eql(u8, opt.type, "mash"); const kmash = std.mem.eql(u8, opt.type, "kmash");
    const lzsh = std.mem.eql(u8, opt.type, "lzsh");

    var pop      = try Matrix(T).init(opt.iterations, mpt.states(opt.potential), allocator); defer      pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations, 1                        , allocator); defer     ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations, 1                        , allocator); defer     epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations, 1                        , allocator); defer     etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer momentum.deinit(); momentum.fill(0);
    var coefs    = try Matrix(T).init(opt.iterations, mpt.states(opt.potential), allocator); defer    coefs.deinit();    coefs.fill(0);

    {
        var T1   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer T1.deinit();
        var T2   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer T2.deinit();

        var r  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  r.deinit();
        var p  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  p.deinit();
        var v  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  v.deinit();
        var a  = try Vector(T).init(mpt.dims(opt.potential), allocator); defer  a.deinit();
        var ap = try Vector(T).init(mpt.dims(opt.potential), allocator); defer ap.deinit();

        var U    = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer   U.deinit();
        var UA   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  UA.deinit();
        var UC   = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer  UC.deinit();
        var UCS  = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer UCS.deinit();
        var TDC  = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer TDC.deinit();

        var C = try Vector(Complex(T)).init(mpt.states(opt.potential), allocator); defer C.deinit(); 
        var S = try Vector(T         ).init(3                        , allocator); defer S.deinit(); 

        var KC1 = try Vector(Complex(T)).init(C.rows, allocator); defer KC1.deinit();
        var KC2 = try Vector(Complex(T)).init(C.rows, allocator); defer KC2.deinit();
        var KC3 = try Vector(Complex(T)).init(C.rows, allocator); defer KC3.deinit();
        var KC4 = try Vector(Complex(T)).init(C.rows, allocator); defer KC4.deinit();
        var KS1 = try Vector(T         ).init(C.rows, allocator); defer KS1.deinit();
        var KS2 = try Vector(T         ).init(C.rows, allocator); defer KS2.deinit();
        var KS3 = try Vector(T         ).init(C.rows, allocator); defer KS3.deinit();
        var KS4 = try Vector(T         ).init(C.rows, allocator); defer KS4.deinit();
        var LS1 = try Vector(T         ).init(C.rows, allocator); defer LS1.deinit();
        var LS2 = try Vector(T         ).init(C.rows, allocator); defer LS2.deinit();
        var LS3 = try Vector(T         ).init(C.rows, allocator); defer LS3.deinit();
        var LS4 = try Vector(T         ).init(C.rows, allocator); defer LS4.deinit();
        var MS1 = try Vector(T         ).init(C.rows, allocator); defer MS1.deinit();
        var MS2 = try Vector(T         ).init(C.rows, allocator); defer MS2.deinit();
        var MS3 = try Vector(T         ).init(C.rows, allocator); defer MS3.deinit();
        var MS4 = try Vector(T         ).init(C.rows, allocator); defer MS4.deinit();

        var U3  = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer  U3[0].deinit(); defer  U3[1].deinit(); defer U3[2].deinit();
        var UC2 = [2]Matrix(T){try U.clone(), try U.clone()               }; defer UC2[0].deinit(); defer UC2[1].deinit()                      ;

        if (print) {try std.io.getStdOut().writer().print("\n{s:6} {s:6} {s:12} {s:12} {s:12} {s:5}", .{"TRAJ", "ITER", "EKIN", "EPOT", "ETOT", "STATE"});}

        if (print) {if (r.rows > 1   )  for (0..r.rows - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}", .{"POSITION"  });}
        if (print) {if (v.rows > 1   )  for (0..v.rows - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}", .{"MOMENTUM"  });}
        if (print) {if (fssh or kfssh) {for (0..C.rows - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}  try std.io.getStdOut().writer().print(" {s:11}", .{"|COEFS|^2" });}}
        if (print) {if (mash or kmash) {for (0..S.rows - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}  try std.io.getStdOut().writer().print(" {s:11}", .{"|BLOCHV|^2"});}}

        if (print) {try std.io.getStdOut().writer().print("\n", .{});}

        for (0..opt.trajectories) |i| {

            const theta = 2 * std.math.pi * rand_bloc.float(T); const phi = if (opt.initial_conditions.state == 1) std.math.acos(1 - rand_bloc.float(T)) else std.math.acos(0 - rand_bloc.float(T));

            for (0..r.rows) |j| r.ptr(j).* = opt.initial_conditions.position_mean[j] + opt.initial_conditions.position_std[j] * rand_traj.floatNorm(T);
            for (0..p.rows) |j| p.ptr(j).* = opt.initial_conditions.momentum_mean[j] + opt.initial_conditions.momentum_std[j] * rand_traj.floatNorm(T);

            if (mash or kmash) {S.ptr(0).* = std.math.sin(phi) * std.math.cos(theta); S.ptr(1).* = std.math.sin(phi) * std.math.sin(theta); S.ptr(2).* = std.math.cos(phi);}

            if (fssh or kfssh) {C.fill(Complex(T).init(0, 0)); C.ptr(opt.initial_conditions.state).* = Complex(T).init(1, 0);}

            @memcpy(v.data, p.data); v.ptr(0).* /= opt.initial_conditions.mass; a.fill(0); var s = opt.initial_conditions.state; var sp = s; var Ekin: T = 0; var Epot: T = 0;

            for (0..opt.iterations) |j| {

                @memcpy(ap.data, a.data); sp = s; 

                for (0..r.rows) |k| {

                    r.ptr(k).* += 1 * opt.derivative_step; mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data);} const Up = U.at(s, s);
                    r.ptr(k).* -= 2 * opt.derivative_step; mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data);} const Um = U.at(s, s);

                    a.ptr(k).* = -0.5 * (Up - Um) / opt.derivative_step / opt.initial_conditions.mass;
                    v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                    r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step + opt.derivative_step;
                }

                mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {

                    mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data); @memcpy(UC2[j % 2].data, UC.data);

                    if (j > 0) for (0..UC.cols) |k| {
                        var overlap: T = 0; for (0..UC.rows) |l| {overlap += UC2[j % 2].at(l, k) * UC2[(j - 1) % 2].at(l, k);} if (overlap < 0) for (0..UC.rows) |l| {UC2[j % 2].ptr(l, k).* *= -1;};
                    };
                }

                @memcpy(U3[j % 3].data, U.data); Ekin = 0; for (v.data) |e| {Ekin += e * e;} Ekin *= 0.5 * opt.initial_conditions.mass; Epot = U.at(s, s);

                if (opt.adiabatic and  (fssh or  mash) and j > 0) derivativeCouplingNumeric(T, &TDC, &UCS, &[_]Matrix(T){UC2[j % 2], UC2[(j - 1) % 2]},                opt.time_step);
                if (opt.adiabatic and (kfssh or kmash) and j > 1) derivativeCouplingBaeckan(T, &TDC,       &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, opt.time_step);

                if ((lzsh         ) and j > 1) s = try landauZener(T, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, opt.adiabatic, rand_jump);
                if ((fssh or kfssh) and j > 1) s = try fewestSwitches(T, &C, U, TDC, s, opt.time_step, rand_jump, &KC1, &KC2, &KC3, &KC4);
                if ((mash or kmash) and j > 1) s = try mappingApproach(T, &S, U, TDC, s, opt.time_step, &KS1, &KS2, &KS3, &KS4, &LS1, &LS2, &LS3, &LS4, &MS1, &MS2, &MS3, &MS4);

                if (s != sp and Ekin < U.at(s, s) - U.at(sp, sp)) s = sp;

                if (sp != s) {for (0..v.rows) |k| v.ptr(k).* *= std.math.sqrt((Ekin - U.at(s, s) + U.at(sp, sp)) / Ekin);}

                if (opt.write.population_mean       != null) pop.ptr(j, s).* += 1;
                if (opt.write.kinetic_energy_mean   != null) ekin.ptr(j, 0).* += Ekin                                                             ;
                if (opt.write.potential_energy_mean != null) epot.ptr(j, 0).* += Epot                                                             ;
                if (opt.write.total_energy_mean     != null) etot.ptr(j, 0).* += Ekin + Epot                                                      ;
                if (opt.write.position_mean         != null) for (0..r.rows) |k| {position.ptr(j, k).* += r.at(k);}                               ;
                if (opt.write.momentum_mean         != null) for (0..v.rows) |k| {momentum.ptr(j, k).* += v.at(k) * opt.initial_conditions.mass;} ;
                if (opt.write.fssh_coefficient_mean != null) for (0..C.rows) |k| {coefs.ptr(j, k).* += C.at(k).magnitude() * C.at(k).magnitude();};

                if (j == opt.iterations - 1) output.pop.ptr(s).* += 1;

                if (print and (i == 0 or (i + 1) % opt.log_intervals.trajectory == 0) and (j == 0 or (j + 1) % opt.log_intervals.iteration == 0)) {
                    try printIteration(T, @intCast(i), @intCast(j), Ekin, Epot, Ekin + Epot, s, r, v, C, S, opt.initial_conditions.mass, fssh or kfssh, mash or kmash);
                }
            }
        }
    }

    for (0..opt.iterations) |i| {for (0..pop.cols     ) |j|     pop.ptr(i, j).*  /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..ekin.cols    ) |j|     ekin.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..epot.cols    ) |j|     epot.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..etot.cols    ) |j|     etot.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..position.cols) |j| position.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..momentum.cols) |j| momentum.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations) |i| {for (0..coefs.cols   ) |j|    coefs.ptr(i, j).* /= asfloat(T, opt.trajectories);}

    for (0..output.pop.rows) |i| output.pop.ptr(i).* /= asfloat(T, opt.trajectories);

    for (0..mpt.states(opt.potential)) |i| {
        if (print) {try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.pop.at(i)});}
    }

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, coefs, allocator); return output;
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

fn mappingApproach(comptime T: type, S: *Vector(T), U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, K1: @TypeOf(S), K2: @TypeOf(S), K3: @TypeOf(S), K4: @TypeOf(S), L1: @TypeOf(S), L2: @TypeOf(S), L3: @TypeOf(S), L4: @TypeOf(S), M1: @TypeOf(S), M2: @TypeOf(S), M3: @TypeOf(S), M4: @TypeOf(S)) !u32 {
    const iters = 10;

    _ = K1; _ = K2; _ = K3; _ = K4;
    _ = L1; _ = L2; _ = L3; _ = L4;
    _ = M1; _ = M2; _ = M3; _ = M4;
    _ = U; _ = TDC; _ = time_step;
    _ = s;

    for (0..iters) |_| {

    }

    return if (S.at(2) > 0) 1 else 0;
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

        if ((di0 - ds0) * (di1 - ds1) > 0) continue;

        const ddi = (di0 - di1) / time_step; const dds = (ds0 - ds1) / time_step;

        const p = std.math.exp(-0.5 * std.math.pi * std.math.sqrt(std.math.pow(T, U3[0].at(i, i) - U3[0].at(s, s), 3) / (ddi - dds)));

        if (!sampled) {rn = rand.float(T); sampled = true;} if (rn < p and p > pm) {sn = @intCast(i); pm = p;}
    };

    return sn;
}

fn printIteration(comptime T: type, i: u32, j: u32, Ekin: T, Epot: T, Etot: T, s: u32, r: Vector(T), v: Vector(T), C: Vector(Complex(T)), S: Vector(T), mass: T, fssh: bool, mash: bool) !void {
    try std.io.getStdOut().writer().print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:5} [", .{i + 1, j + 1, Ekin, Epot, Etot, s});

    for (0..r.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", r.at(k)});
    }

    try std.io.getStdOut().writer().print("] [", .{});

    for (0..v.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", v.at(k) * mass});
    }

    if (fssh or mash) try std.io.getStdOut().writer().print("] [", .{});

    if (fssh) for (0..C.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", C.at(k).magnitude() * C.at(k).magnitude()});
    };

    if (mash) for (0..S.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", S.at(k)});
    };

    try std.io.getStdOut().writer().print("]\n", .{});
}

fn writeResults(comptime T: type, opt: ClassicalDynamicsOptions(T), pop: Matrix(T), ekin: Matrix(T), epot: Matrix(T), etot: Matrix(T), position: Matrix(T), momentum: Matrix(T), coefs: Matrix(T), allocator: std.mem.Allocator) !void {
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

    if (opt.write.fssh_coefficient_mean) |path| {
        var coefs_t = try Matrix(T).init(opt.iterations, coefs.cols + 1, allocator); mat.hjoin(T, &coefs_t, time, coefs); try coefs_t.write(path); coefs_t.deinit();
    }
}

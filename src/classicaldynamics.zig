//! Module for classical dynamics simulations.

const std = @import("std"); const Complex = std.math.Complex;

const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;
const     max = @import("helper.zig").max    ;
const     min = @import("helper.zig").min    ;

/// The classical dynamics options.
pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            position_mean: []const T,
            position_std:  []const T,
            momentum_mean: []const T,
            momentum_std:  []const T,
            state: u32, mass: T
        };
        pub const FewestSwitches = struct {
            quantum_substep: u32 = 10,
            decoherence_alpha: ?T = null
        };
        pub const LandauZener = struct {
        };
        pub const SpinMapping = struct {
            quantum_substep: u32 = 10,
            sample_bloch_vector: bool = true
        };
        pub const LogIntervals = struct {
            trajectory: u32 = 1, iteration: u32 = 1
        };
        pub const Write = struct {
            fssh_coefficient_mean:         ?[]const u8 = null,
            kinetic_energy_mean:           ?[]const u8 = null,
            momentum_mean:                 ?[]const u8 = null,
            population_mean:               ?[]const u8 = null,
            position_mean:                 ?[]const u8 = null,
            potential_energy_mean:         ?[]const u8 = null,
            time_derivative_coupling_mean: ?[]const u8 = null,
            total_energy_mean:             ?[]const u8 = null
        };

        adiabatic: bool,
        derivative_step: T = 0.001,
        iterations: u32,
        potential: []const u8,
        seed: u32 = 1,
        time_derivative_coupling: ?[]const u8 = "numeric",
        time_step: T,
        trajectories: u32,

        fewest_switches: ?FewestSwitches = null, landau_zener: ?LandauZener = null, spin_mapping: ?SpinMapping = null,

        initial_conditions: InitialConditions, log_intervals: LogIntervals = .{}, write: Write = .{},
    };
}

/// The classical dynamics output.
pub fn ClassicalDynamicsOutput(comptime T: type) type {
    return struct {
        pop: Vector(T),

        r: Vector(T),
        p: Vector(T),

        Ekin: T,
        Epot: T,

        /// Initialize the classical dynamics output.
        pub fn init(ndim: usize, nstate: usize, allocator: std.mem.Allocator) !ClassicalDynamicsOutput(T) {
            return ClassicalDynamicsOutput(T){
                .pop = try Vector(T).init(nstate, allocator),

                .r = try Vector(T).init(ndim, allocator),
                .p = try Vector(T).init(ndim, allocator),

                .Ekin = 0,
                .Epot = 0
            };
        }

        /// Free the memory allocated for the classical dynamics output.
        pub fn deinit(self: ClassicalDynamicsOutput(T)) void {
            self.pop.deinit(); self.r.deinit(); self.p.deinit();
        }
    };
}

/// Main function for running classical dynamics simulations.
pub fn run(comptime T: type, opt: ClassicalDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !ClassicalDynamicsOutput(T) {
    if (opt.initial_conditions.state > try mpt.states(opt.potential) - 1) return error.InvalidInitialState;

    const ndim = try mpt.dims(opt.potential); const nstate = try mpt.states(opt.potential);

    var output = try ClassicalDynamicsOutput(T).init(ndim, nstate, allocator);

    var prng_jump = std.Random.DefaultPrng.init(opt.seed); const rand_jump = prng_jump.random();
    var prng_traj = std.Random.DefaultPrng.init(opt.seed); const rand_traj = prng_traj.random();
    var prng_bloc = std.Random.DefaultPrng.init(opt.seed); const rand_bloc = prng_bloc.random();

    var pop      = try Matrix(T).init(opt.iterations + 1, 1 + nstate         , allocator); defer      pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations + 1, 1 + 1              , allocator); defer     ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations + 1, 1 + 1              , allocator); defer     epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations + 1, 1 + 1              , allocator); defer     etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations + 1, 1 + ndim           , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations + 1, 1 + ndim           , allocator); defer momentum.deinit(); momentum.fill(0);
    var coef     = try Matrix(T).init(opt.iterations + 1, 1 + nstate         , allocator); defer     coef.deinit();     coef.fill(0);
    var tdc      = try Matrix(T).init(opt.iterations + 1, 1 + nstate * nstate, allocator); defer      tdc.deinit();      tdc.fill(0);

    {
        const time = try Matrix(T).init(opt.iterations + 1, 1, allocator); defer time.deinit(); time.linspace(0, opt.time_step * asfloat(T, opt.iterations));

        for (0..opt.iterations + 1) |i| {
            pop.ptr(i, 0).*      = time.at(i, 0);
            ekin.ptr(i, 0).*     = time.at(i, 0);
            epot.ptr(i, 0).*     = time.at(i, 0);
            etot.ptr(i, 0).*     = time.at(i, 0);
            position.ptr(i, 0).* = time.at(i, 0);
            momentum.ptr(i, 0).* = time.at(i, 0);
            coef.ptr(i, 0).*     = time.at(i, 0);
            tdc.ptr(i, 0).*      = time.at(i, 0);
        }
    }

    const fssh = opt.fewest_switches != null;
    const lzsh = opt.landau_zener    != null;
    const mash = opt.spin_mapping    != null;

    if ((fssh and lzsh and mash) or (fssh and lzsh) or (fssh and mash) or (lzsh and mash)) {
        return error.MultipleHoppingMechanisms;
    }

    const tdc_numeric = std.mem.eql(u8, opt.time_derivative_coupling.?, "numeric");
    const tdc_baeckan = std.mem.eql(u8, opt.time_derivative_coupling.?, "baeckan");

    if (opt.time_derivative_coupling != null and !tdc_numeric and !tdc_baeckan) {
        return error.UnknownTimeDerivativeCoupling;
    }

    {
        var T1 = try Matrix(T).init(nstate, nstate, allocator); defer T1.deinit();
        var T2 = try Matrix(T).init(nstate, nstate, allocator); defer T2.deinit();

        var r  = try Vector(T).init(ndim, allocator); defer  r.deinit();
        var p  = try Vector(T).init(ndim, allocator); defer  p.deinit();
        var v  = try Vector(T).init(ndim, allocator); defer  v.deinit();
        var a  = try Vector(T).init(ndim, allocator); defer  a.deinit();
        var ap = try Vector(T).init(ndim, allocator); defer ap.deinit();

        var U   = try Matrix(T).init(nstate, nstate, allocator); defer   U.deinit();
        var UA  = try Matrix(T).init(nstate, nstate, allocator); defer  UA.deinit();
        var UC  = try Matrix(T).init(nstate, nstate, allocator); defer  UC.deinit();
        var UCS = try Matrix(T).init(nstate, nstate, allocator); defer UCS.deinit();
        var TDC = try Matrix(T).init(nstate, nstate, allocator); defer TDC.deinit();

        var P = try Vector(T         ).init(nstate, allocator); defer P.deinit();
        var C = try Vector(Complex(T)).init(nstate, allocator); defer C.deinit(); 
        var S = try Vector(T         ).init(3     , allocator); defer S.deinit(); 

        var KC1 = try Vector(Complex(T)).init(C.rows, allocator); defer KC1.deinit();
        var KC2 = try Vector(Complex(T)).init(C.rows, allocator); defer KC2.deinit();
        var KC3 = try Vector(Complex(T)).init(C.rows, allocator); defer KC3.deinit();
        var KC4 = try Vector(Complex(T)).init(C.rows, allocator); defer KC4.deinit();

        var U3  = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer  U3[0].deinit(); defer  U3[1].deinit(); defer U3[2].deinit();
        var UC2 = [2]Matrix(T){try U.clone(), try U.clone()               }; defer UC2[0].deinit(); defer UC2[1].deinit()                      ;

        if (print) {try std.io.getStdOut().writer().print("\n{s:6} {s:6} {s:12} {s:12} {s:12} {s:5}", .{"TRAJ", "ITER", "EKIN", "EPOT", "ETOT", "STATE"});}

        if (print and r.rows > 1) for (0..min(usize, r.rows, 3) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});};

        if (print and r.rows > 3) try std.io.getStdOut().writer().print("     ", .{}); if (print) try std.io.getStdOut().writer().print(" {s:11}", .{"POSITION"  });

        if (print and v.rows > 1) for (0..min(usize, v.rows, 3) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});};

        if (print and v.rows > 3) try std.io.getStdOut().writer().print("     ", .{}); if (print) try std.io.getStdOut().writer().print(" {s:11}", .{"MOMENTUM"  });

        if (print) {if (fssh      ) {for (0..min(usize, C.rows, 4) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}  try std.io.getStdOut().writer().print(" {s:11}", .{"|COEFS|^2" });}}
        if (print) {if (mash      ) {for (0..S.rows - 1               ) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}  try std.io.getStdOut().writer().print(" {s:11}", .{"|BLOCHV|^2"});}}

        if (print) {try std.io.getStdOut().writer().print("\n", .{});}

        for (0..opt.trajectories) |i| {

            for (0..r.rows) |j| r.ptr(j).* = opt.initial_conditions.position_mean[j] + opt.initial_conditions.position_std[j] * rand_traj.floatNorm(T);
            for (0..p.rows) |j| p.ptr(j).* = opt.initial_conditions.momentum_mean[j] + opt.initial_conditions.momentum_std[j] * rand_traj.floatNorm(T);

            if (fssh) {C.fill(Complex(T).init(0, 0)); C.ptr(opt.initial_conditions.state).* = Complex(T).init(1, 0);         }
            if (mash) {try initialBlochVector(T, &S, opt.spin_mapping.?, opt.initial_conditions.state, rand_bloc, allocator);}

            @memcpy(v.data, p.data); for (v.data) |*e| e.* /= opt.initial_conditions.mass; a.fill(0); var s = opt.initial_conditions.state; var sp = s; var Ekin: T = 0; var Epot: T = 0;

            for (0..opt.iterations + 1) |j| {

                @memcpy(ap.data, a.data); sp = s; 

                if (j > 0) for (0..r.rows) |k| {

                    r.ptr(k).* += 1 * opt.derivative_step; try mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data);} const Up = U.at(s, s);
                    r.ptr(k).* -= 2 * opt.derivative_step; try mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data);} const Um = U.at(s, s);

                    a.ptr(k).* = -0.5 * (Up - Um) / opt.derivative_step / opt.initial_conditions.mass;
                    v.ptr(k).* += 0.5 * (a.at(k) + ap.at(k)) * opt.time_step;
                    r.ptr(k).* += (v.at(k) + 0.5 * a.at(k) * opt.time_step) * opt.time_step + opt.derivative_step;
                };

                try mpt.eval(T, &U, opt.potential, r); if (opt.adiabatic) {

                    mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data); @memcpy(UC2[j % 2].data, UC.data);

                    if (j > 0) for (0..UC.cols) |k| {
                        var overlap: T = 0; for (0..UC.rows) |l| {overlap += UC2[j % 2].at(l, k) * UC2[(j - 1) % 2].at(l, k);} if (overlap < 0) for (0..UC.rows) |l| {UC2[j % 2].ptr(l, k).* *= -1;};
                    };
                }

                @memcpy(U3[j % 3].data, U.data); Ekin = 0; for (v.data) |e| Ekin += e * e; Ekin *= 0.5 * opt.initial_conditions.mass; Epot = U.at(s, s);

                if (opt.adiabatic and tdc_numeric and j > 0) derivativeCouplingNumeric(T, &TDC, &UCS, &[_]Matrix(T){UC2[j % 2], UC2[(j - 1) % 2]},                opt.time_step);
                if (opt.adiabatic and tdc_baeckan and j > 1) derivativeCouplingBaeckan(T, &TDC,       &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, opt.time_step);

                if (lzsh and j > 1) s = landauZener(T, &P, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, opt.adiabatic, rand_jump);
                if (fssh and j > 1) s = try fewestSwitches(T, &C, opt.fewest_switches.?, U, TDC, s, opt.time_step, Ekin, rand_jump, &KC1, &KC2, &KC3, &KC4);
                if (mash and j > 1) s = try spinMapping(T, &S, opt.spin_mapping.?, U, TDC, s, opt.time_step, rand_jump);

                if (s != sp and Ekin < U.at(s, s) - U.at(sp, sp)) s = sp;

                if (sp != s) for (0..v.rows) |k| {
                    v.ptr(k).* *= std.math.sqrt((Ekin - U.at(s, s) + U.at(sp, sp)) / Ekin);
                };

                if (opt.write.population_mean               != null) pop.ptr(j, 1 + s).* += 1;
                if (opt.write.kinetic_energy_mean           != null) ekin.ptr(j, 1 + 0).* += Ekin                                                            ;
                if (opt.write.potential_energy_mean         != null) epot.ptr(j, 1 + 0).* += Epot                                                            ;
                if (opt.write.total_energy_mean             != null) etot.ptr(j, 1 + 0).* += Ekin + Epot                                                     ;
                if (opt.write.position_mean                 != null) for (0..r.rows) |k| {position.ptr(j, 1 + k).* += r.at(k);}                              ;
                if (opt.write.momentum_mean                 != null) for (0..v.rows) |k| {momentum.ptr(j, 1 + k).* += v.at(k) * opt.initial_conditions.mass;};
                if (opt.write.time_derivative_coupling_mean != null) for (0..TDC.rows * TDC.cols) |k| {tdc.ptr(j, 1 + k).* += TDC.data[k];}                  ;
                if (opt.write.fssh_coefficient_mean         != null) for (0..C.rows) |k| {coef.ptr(j, 1 + k).* += C.at(k).magnitude() * C.at(k).magnitude();};

                if (j == opt.iterations - 1) {

                    output.pop.ptr(s).* += 1;

                    for (0..ndim) |k| {
                        output.r.ptr(k).* += r.at(k);  output.p.ptr(k).* += v.at(k) * opt.initial_conditions.mass;
                    }

                    output.Ekin += Ekin;
                    output.Epot += Epot;
                }

                if (print and (i == 0 or (i + 1) % opt.log_intervals.trajectory == 0) and (j % opt.log_intervals.iteration == 0)) {
                    try printIteration(T, @intCast(i), @intCast(j), Ekin, Epot, Ekin + Epot, s, r, v, C, S, opt.initial_conditions.mass, fssh, mash);
                }
            }
        }
    }

    for (0..opt.iterations + 1) |i| {for (1..pop.cols     ) |j|      pop.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..ekin.cols    ) |j|     ekin.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..epot.cols    ) |j|     epot.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..etot.cols    ) |j|     etot.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..position.cols) |j| position.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..momentum.cols) |j| momentum.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..tdc.cols     ) |j|      tdc.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..coef.cols    ) |j|     coef.ptr(i, j).* /= asfloat(T, opt.trajectories);}

    for (0..output.pop.rows) |i| output.pop.ptr(i).* /= asfloat(T, opt.trajectories);
    for (0..output.r.rows  ) |i|   output.r.ptr(i).* /= asfloat(T, opt.trajectories);
    for (0..output.p.rows  ) |i|   output.p.ptr(i).* /= asfloat(T, opt.trajectories);

    output.Ekin /= asfloat(T, opt.trajectories);
    output.Epot /= asfloat(T, opt.trajectories);

    for (0..nstate) |i| {
        if (print) {try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.pop.at(i)});}
    }

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, tdc, coef); return output;
}

/// Calculate the nonadiabatic coupling between two states using the Baeck-An method.
pub fn derivativeCouplingBaeckan(comptime T: type, TDC: *Matrix(T), U3: []const Matrix(T), time_step: T) void {
    TDC.fill(0);

    for (0..TDC.rows) |i| for (i + 1..TDC.cols) |j| {

        const di0 = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const dj0 = (U3[0].at(j, j) - U3[1].at(j, j)) / time_step;
        const di1 = (U3[1].at(i, i) - U3[2].at(i, i)) / time_step; const dj1 = (U3[1].at(j, j) - U3[2].at(j, j)) / time_step;

        const ddi = (di0 - di1) / time_step; const ddj = (dj0 - dj1) / time_step; const arg = (ddi - ddj) / (U3[0].at(i, i) - U3[0].at(j, j));

        if (arg > 0) {TDC.ptr(i, j).* = 0.5 * std.math.sqrt(arg);} TDC.ptr(j, i).* = -TDC.at(i, j);
    };
}

/// Calculate the nonadiabatic coupling between two states numerically.
pub fn derivativeCouplingNumeric(comptime T: type, TDC: *Matrix(T), UCS: *Matrix(T), UC2: []const Matrix(T), time_step: T) void {
    UCS.fill(0);

    for (0..UCS.rows) |i| for (0..UCS.cols) |j| for (0..UCS.rows) |k| {
        UCS.ptr(i, j).* += UC2[1].at(k, i) * UC2[0].at(k, j);
    };

    for (0..TDC.rows) |i| for (0..TDC.cols) |j| {
        TDC.ptr(i, j).* = (UCS.at(i, j) - UCS.at(j, i)) / (2 * time_step);
    };
}

/// Function to propagate the wavefunction coefficients used in the FSSH method. The function returns the new state, if a switch occurs.
pub fn fewestSwitches(comptime T: type, C: *Vector(Complex(T)), opt: ClassicalDynamicsOptions(T).FewestSwitches, U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, Ekin: T, rand: std.Random, K1: @TypeOf(C), K2: @TypeOf(C), K3: @TypeOf(C), K4: @TypeOf(C)) !u32 {
    const iters = asfloat(T, opt.quantum_substep); const alpha = if (opt.decoherence_alpha == null) std.math.inf(T) else opt.decoherence_alpha.?; var ns = s;

    const Function = struct { fn get (K: *Vector(Complex(T)), FC: Vector(Complex(T)), FU: Matrix(T), FTDC: Matrix(T)) void {
        for (0..FC.rows) |i| {K.ptr(i).* = FC.at(i).mul(Complex(T).init(FU.at(i, i), 0)).mulbyi().neg();} for (0..FC.rows) |i| for (0..FC.rows) |j| {
            K.ptr(i).* = K.at(i).sub(FC.at(j).mul(Complex(T).init(FTDC.at(i, j), 0)));
        };
    }};

    for (0..opt.quantum_substep) |_| {

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

        const rn = rand.float(T); var p: T = 0; for (0..C.rows) |j| if (j != ns) {
            p += 2 * TDC.at(ns, j) * C.at(j).mul(C.at(ns).conjugate()).re / std.math.pow(T, C.at(ns).magnitude(), 2) * time_step / iters; if (rn < p) {ns = @intCast(j); break;}
        };

        for (0..C.rows) |j| if (j != ns) {
            C.ptr(j).* = C.at(j).mul(Complex(T).init(std.math.exp(-0.5 * @abs(U.at(j, j) - U.at(ns, ns)) * time_step / iters / (1 + alpha / Ekin)), 0));
        };

        var sumc: T = 0; for (0..C.rows) |j| if (j != ns) {sumc += C.at(j).magnitude() * C.at(j).magnitude();};

        if (C.at(ns).magnitude() > 0) {
            C.ptr(ns).* = C.at(ns).mul(Complex(T).init(std.math.sqrt((1 - sumc) / C.at(ns).magnitude() / C.at(ns).magnitude()), 0));
        }
    }

    return ns;
}

/// Function to initialize the initial vector on the Bloch sphere for the spin mapping methods.
pub fn initialBlochVector(comptime T: type, S: *Vector(T), opt: ClassicalDynamicsOptions(T).SpinMapping, s: u32, rand: std.Random, allocator: std.mem.Allocator) !void {
    const theta: T = if (s == 1) 0 else std.math.pi; const phi: T = 0; const xi = 2 * std.math.pi * rand.float(T);

    S.ptr(0).* = std.math.sin(theta) * std.math.cos(phi); S.ptr(1).* = std.math.sin(theta) * std.math.sin(phi); S.ptr(2).* = std.math.cos(theta);

    if (opt.sample_bloch_vector) {

        const x = try Vector(T).init(3, allocator); defer x.deinit();
        const y = try Vector(T).init(3, allocator); defer y.deinit();

        x.ptr(0).* = std.math.cos(theta) * std.math.cos(phi); x.ptr(1).* = std.math.cos(theta) * std.math.sin(phi); x.ptr(2).* = -std.math.sin(theta);
        y.ptr(0).* = -std.math.sin(phi);                      y.ptr(1).* = std.math.cos(phi);                       y.ptr(2).* = 0;

        for (0..S.rows) |j| {
            S.ptr(j).* = S.at(j) / std.math.sqrt(3) + (x.at(j) * std.math.cos(xi) + y.at(j) * std.math.sin(xi)) * std.math.sqrt2 / std.math.sqrt(3);
        }
    }
}

/// Function to propagate the vector on the Bloch sphere for the spin mapping methods. The function returns the new state, if a switch occurs.
pub fn spinMapping(comptime T: type, S: *Vector(T), opt: ClassicalDynamicsOptions(T).SpinMapping, U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, rand: std.Random) !u32 {
    const iters = asfloat(T, opt.quantum_substep); var sn = s;

    const FunctionX = struct { fn get (FS: Vector(T), FU: Matrix(T), FTDC: Matrix(T)) T {
        return (FU.at(1, 1) - FU.at(0, 0)) * FS.at(1) - 2 * FTDC.at(0, 1) * FS.at(2);
    }};

    const FunctionY = struct { fn get (FS: Vector(T), FU: Matrix(T)                 ) T {
        return -(FU.at(1, 1) - FU.at(0, 0)) * FS.at(0);
    }};

    const FunctionZ = struct { fn get (FS: Vector(T),                FTDC: Matrix(T)) T {
        return 2 * FTDC.at(0, 1) * FS.at(0);
    }};

    for (0..opt.quantum_substep) |_| {

        const kx1 = FunctionX.get(S.*, U, TDC);
        const ky1 = FunctionY.get(S.*, U,    );
        const kz1 = FunctionZ.get(S.*,    TDC);

        (S.ptr(0)).* += 0.5 * kx1 * time_step / iters;
        (S.ptr(1)).* += 0.5 * ky1 * time_step / iters;
        (S.ptr(2)).* += 0.5 * kz1 * time_step / iters;

        const kx2 = FunctionX.get(S.*, U, TDC);
        const ky2 = FunctionY.get(S.*, U     );
        const kz2 = FunctionZ.get(S.*,    TDC);

        S.ptr(0).* += 0.5 * (kx2 - kx1) * time_step / iters;
        S.ptr(1).* += 0.5 * (ky2 - ky1) * time_step / iters;
        S.ptr(2).* += 0.5 * (kz2 - kz1) * time_step / iters;

        const kx3 = FunctionX.get(S.*, U, TDC);
        const ky3 = FunctionY.get(S.*, U     );
        const kz3 = FunctionZ.get(S.*,    TDC);

        S.ptr(0).* += 0.5 * (2 * kx3 - kx2) * time_step / iters;
        S.ptr(1).* += 0.5 * (2 * ky3 - ky2) * time_step / iters;
        S.ptr(2).* += 0.5 * (2 * kz3 - kz2) * time_step / iters;

        const kx4 = FunctionX.get(S.*, U, TDC);
        const ky4 = FunctionY.get(S.*, U     );
        const kz4 = FunctionZ.get(S.*,    TDC);

        S.ptr(0).* -= kx3 * time_step / iters;
        S.ptr(1).* -= ky3 * time_step / iters;
        S.ptr(2).* -= kz3 * time_step / iters;

        S.ptr(0).* += time_step / iters / 6 * (kx1 + 2 * kx2 + 2 * kx3 + kx4);
        S.ptr(1).* += time_step / iters / 6 * (ky1 + 2 * ky2 + 2 * ky3 + ky4);
        S.ptr(2).* += time_step / iters / 6 * (kz1 + 2 * kz2 + 2 * kz3 + kz4);

        const rn = rand.float(T); if (!opt.sample_bloch_vector) {
            if (s == sn and s == 0 and rn <  FunctionZ.get(S.*, TDC) / (1 - S.at(2)) * time_step / iters) sn = 1;
            if (s == sn and s == 1 and rn < -FunctionZ.get(S.*, TDC) / (1 + S.at(2)) * time_step / iters) sn = 0;
        }
    }

    if (opt.sample_bloch_vector) {return if (S.at(2) > 0) 1 else 0;} else return sn;
}

/// Function to calculate the Landau-Zener probability of a transition between two states. The function returns the new state, if a switch occurs.
pub fn landauZener(comptime T: type, P: *Vector(T), U3: []const Matrix(T), s: u32, time_step: T, adiabatic: bool, rand: std.Random) u32 {
    var sn = s; var rn: T = 0; P.fill(0);

    if (!adiabatic) for (0..U3[0].rows) |i| if (i != s) {

        if ((U3[0].at(i, i) - U3[0].at(s, s)) * (U3[1].at(i, i) - U3[1].at(s, s)) > 0) continue;

        const di = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const ds = (U3[0].at(s, s) - U3[1].at(s, s)) / time_step;

        P.ptr(i).* = 1 - std.math.exp(-2 * std.math.pi * std.math.pow(T, U3[0].at(i, s), 2) / @abs(di - ds));

        if (std.math.isNan(P.at(i))) P.ptr(i).* = 0;
    };

    if (adiabatic) for (0..U3[0].rows) |i| if (i != s) {

        const di0 = (U3[0].at(i, i) - U3[1].at(i, i)) / time_step; const ds0 = (U3[0].at(s, s) - U3[1].at(s, s)) / time_step;
        const di1 = (U3[1].at(i, i) - U3[2].at(i, i)) / time_step; const ds1 = (U3[1].at(s, s) - U3[2].at(s, s)) / time_step;

        if ((di0 - ds0) * (di1 - ds1) > 0) continue;

        const ddi = (di0 - di1) / time_step; const dds = (ds0 - ds1) / time_step;

        P.ptr(i).* = std.math.exp(-0.5 * std.math.pi * std.math.sqrt(std.math.pow(T, U3[0].at(i, i) - U3[0].at(s, s), 3) / (ddi - dds)));

        if (std.math.isNan(P.at(i))) P.ptr(i).* = 0;
    };

    var ps: T = 0; for (0..P.rows) |i| ps += P.at(i);

    if (ps > 0) {

        if (ps > 1) {
            for (0..P.rows) |i| P.ptr(i).* /= ps;
        }

        for (1..P.rows) |i| P.ptr(i).* += P.at(i - 1);

        rn = rand.float(T);

        for (0..P.rows) |i| if (rn > (if (i > 0) P.at(i - 1) else 0) and rn < P.at(i)) {
            sn = @intCast(i); break;
        };
    }

    return sn;
}

/// Function to print the results of a single iteration.
pub fn printIteration(comptime T: type, i: u32, j: u32, Ekin: T, Epot: T, Etot: T, s: u32, r: Vector(T), v: Vector(T), C: Vector(Complex(T)), S: Vector(T), mass: T, fssh: bool, mash: bool) !void {
    try std.io.getStdOut().writer().print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:5} [", .{i + 1, j, Ekin, Epot, Etot, s});

    for (0..min(usize, r.rows, 3)) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", r.at(k)});
    }

    if (r.rows > 3) {try std.io.getStdOut().writer().print(", ...", .{});}

    try std.io.getStdOut().writer().print("] [", .{});

    for (0..min(usize, v.rows, 3)) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", v.at(k) * mass});
    }

    if (v.rows > 3) {try std.io.getStdOut().writer().print(", ...", .{});}

    if (fssh or mash) try std.io.getStdOut().writer().print("] [", .{});

    if (fssh) for (0..C.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", C.at(k).magnitude() * C.at(k).magnitude()});
    };

    if (mash) for (0..S.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", S.at(k)});
    };

    try std.io.getStdOut().writer().print("]\n", .{});
}

/// Function to write the results to disk.
pub fn writeResults(comptime T: type, opt: ClassicalDynamicsOptions(T), pop: Matrix(T), ekin: Matrix(T), epot: Matrix(T), etot: Matrix(T), position: Matrix(T), momentum: Matrix(T), tdc: Matrix(T), coef: Matrix(T)) !void {
    if (opt.write.kinetic_energy_mean          ) |path| try     ekin.write(path);
    if (opt.write.population_mean              ) |path| try      pop.write(path);
    if (opt.write.potential_energy_mean        ) |path| try     epot.write(path);
    if (opt.write.total_energy_mean            ) |path| try     etot.write(path);
    if (opt.write.position_mean                ) |path| try position.write(path);
    if (opt.write.momentum_mean                ) |path| try momentum.write(path);
    if (opt.write.time_derivative_coupling_mean) |path| try      tdc.write(path);
    if (opt.write.fssh_coefficient_mean        ) |path| try      coef.write(path);
}

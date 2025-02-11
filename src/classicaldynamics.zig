//! Module for classical dynamics simulations.

const std = @import("std"); const Complex = std.math.Complex;

const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");
const vec = @import("vector.zig"        );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const contains = @import("helper.zig").contains;
const      abs = @import("helper.zig").abs     ;
const      max = @import("helper.zig").max     ;
const      min = @import("helper.zig").min     ;
const      sgn = @import("helper.zig").sgn     ;
const      sum = @import("helper.zig").sum     ;
const  asfloat = @import("helper.zig").asfloat ;

/// The classical dynamics options.
pub fn ClassicalDynamicsOptions(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            position_mean: []const T,
            position_std:  []const T,
            momentum_mean: []const T,
            momentum_std:  []const T,
            state:         []const T,
            mass:          []const T
        };
        pub const FewestSwitches = struct {
            quantum_substep: u32 = 10,
            decoherence_alpha: ?T = null
        };
        pub const LandauZener = struct {
        };
        pub const SpinMapping = struct {
            fewest_switches: bool = false, quantum_jump_iteration: ?[]const u32 = null
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
    const ndim = try mpt.dims(opt.potential); const nstate = try mpt.states(opt.potential);

    if (opt.initial_conditions.state.len         != nstate) return error.InvalidInitialState;
    if (opt.initial_conditions.position_mean.len != ndim  ) return error.InvalidMeanPosition;
    if (opt.initial_conditions.position_std.len  != ndim  ) return error.InvalidStdPosition ;
    if (opt.initial_conditions.momentum_mean.len != ndim  ) return error.InvalidMeanMomentum;
    if (opt.initial_conditions.momentum_std.len  != ndim  ) return error.InvalidStdMomentum ;
    if (opt.initial_conditions.mass.len          != ndim  ) return error.InvalidMass        ;
    if (sum(T, opt.initial_conditions.state)     != 1     ) return error.InvalidStateSum    ;

    var output = try ClassicalDynamicsOutput(T).init(ndim, nstate, allocator);

    var prng_jump = std.Random.DefaultPrng.init(opt.seed); const rand_jump = prng_jump.random();
    var prng_traj = std.Random.DefaultPrng.init(opt.seed); const rand_traj = prng_traj.random();
    var prng_bloc = std.Random.DefaultPrng.init(opt.seed); const rand_bloc = prng_bloc.random();

    var coef     = try Matrix(T).init(opt.iterations + 1, 1 + nstate         , allocator); defer     coef.deinit(); coef    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var ekin     = try Matrix(T).init(opt.iterations + 1, 1 + 1              , allocator); defer     ekin.deinit(); ekin    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var epot     = try Matrix(T).init(opt.iterations + 1, 1 + 1              , allocator); defer     epot.deinit(); epot    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var etot     = try Matrix(T).init(opt.iterations + 1, 1 + 1              , allocator); defer     etot.deinit(); etot    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var momentum = try Matrix(T).init(opt.iterations + 1, 1 + ndim           , allocator); defer momentum.deinit(); momentum.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var pop      = try Matrix(T).init(opt.iterations + 1, 1 + nstate         , allocator); defer      pop.deinit(); pop     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var position = try Matrix(T).init(opt.iterations + 1, 1 + ndim           , allocator); defer position.deinit(); position.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var tdc      = try Matrix(T).init(opt.iterations + 1, 1 + nstate * nstate, allocator); defer      tdc.deinit(); tdc     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));

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

        var r = try Vector(T).init(ndim, allocator); defer r.deinit();
        var p = try Vector(T).init(ndim, allocator); defer p.deinit();
        var v = try Vector(T).init(ndim, allocator); defer v.deinit();
        var a = try Vector(T).init(ndim, allocator); defer a.deinit();

        var U   = try Matrix(T).init(nstate, nstate, allocator); defer   U.deinit();
        var UA  = try Matrix(T).init(nstate, nstate, allocator); defer  UA.deinit();
        var UC  = try Matrix(T).init(nstate, nstate, allocator); defer  UC.deinit();
        var UCS = try Matrix(T).init(nstate, nstate, allocator); defer UCS.deinit();
        var TDC = try Matrix(T).init(nstate, nstate, allocator); defer TDC.deinit();

        var P  = try Vector(T         ).init(nstate    ,    allocator); defer P.deinit();
        var C  = try Vector(Complex(T)).init(nstate    ,    allocator); defer C.deinit();
        var S  = try Matrix(T         ).init(nstate - 1, 3, allocator); defer S.deinit();
        var SI = try Matrix(usize     ).init(nstate - 1, 2, allocator); defer SI.deinit();

        var KC1 = try Vector(Complex(T)).init(C.rows, allocator); defer KC1.deinit();
        var KC2 = try Vector(Complex(T)).init(C.rows, allocator); defer KC2.deinit();
        var KC3 = try Vector(Complex(T)).init(C.rows, allocator); defer KC3.deinit();
        var KC4 = try Vector(Complex(T)).init(C.rows, allocator); defer KC4.deinit();

        var SP = try Matrix(T).init(3,          3, allocator); defer SP.deinit();
        var SN = try Matrix(T).init(nstate - 1, 3, allocator); defer SN.deinit();
        var S0 = try Matrix(T).init(nstate - 1, 3, allocator); defer S0.deinit();

        var U3  = [3]Matrix(T){try U.clone(), try U.clone(), try U.clone()}; defer  U3[0].deinit(); defer  U3[1].deinit(); defer U3[2].deinit();
        var UC2 = [2]Matrix(T){try U.clone(), try U.clone()               }; defer UC2[0].deinit(); defer UC2[1].deinit()                      ;

        if (print) try printHeader(ndim, nstate, fssh, mash);

        for (0..opt.trajectories) |i| {

            var s: u32 = undefined; for (0..nstate) |j| if (asfloat(T, i + 1) / asfloat(T, opt.trajectories) <= sum(T, opt.initial_conditions.state[0..j + 1])) {s = @intCast(j); break;};

            for (0..r.rows) |j| r.ptr(j).* = opt.initial_conditions.position_mean[j] + opt.initial_conditions.position_std[j] * rand_traj.floatNorm(T);
            for (0..p.rows) |j| p.ptr(j).* = opt.initial_conditions.momentum_mean[j] + opt.initial_conditions.momentum_std[j] * rand_traj.floatNorm(T);

            a.fill(0); @memcpy(v.data, p.data); for (0..v.rows) |j| v.ptr(j).* /= opt.initial_conditions.mass[j];

            if (fssh) {C.fill(Complex(T).init(0, 0)); C.ptr(s).* = Complex(T).init(1, 0);}
            if (mash) {try initialBlochVector(T, &S, &SI, opt.spin_mapping.?, s, rand_bloc);  }

            @memcpy(S0.data, S.data); var ns = s; var Ekin: T = 0; var Epot: T = 0;

            for (0..opt.iterations + 1) |j| {

                if (j > 0) try propagate(T, opt, &r, &v, &a, &U, &UA, &UC, &T1, &T2, s);

                try mpt.eval(T, &U, opt.potential, r);

                if (opt.adiabatic) adiabatizePotential(T, &U, &UA, &UC, &UC2, &T1, &T2, j);

                @memcpy(U3[j % 3].data, U.data); Ekin = 0; for (0..v.rows) |k| Ekin += 0.5 * opt.initial_conditions.mass[k] * v.at(k) * v.at(k); Epot = U.at(s, s);

                if (opt.adiabatic and tdc_numeric and j > 0) derivativeCouplingNumeric(T, &TDC, &UCS, &[_]Matrix(T){UC2[j % 2], UC2[(j - 1) % 2]},                opt.time_step);
                if (opt.adiabatic and tdc_baeckan and j > 1) derivativeCouplingBaeckan(T, &TDC,       &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, opt.time_step);

                if (lzsh and j > 1) ns = landauZener(T, &P, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, opt.adiabatic, rand_jump);
                if (fssh and j > 1) ns = try fewestSwitches(T, &C, opt.fewest_switches.?, U, TDC, s, opt.time_step, Ekin, rand_jump, &KC1, &KC2, &KC3, &KC4);
                if (mash and j > 1) ns = try spinMapping(T, &S, &SI, opt.spin_mapping.?, U, TDC, s, opt.time_step, rand_jump, &SP, &SN);

                if (s != ns and Ekin >= U.at(ns, ns) - U.at(s, s)) {
                    rescaleVelocity(T, &v, s, ns, U, Ekin); s = ns;
                }

                if (mash and opt.spin_mapping.?.quantum_jump_iteration != null and contains(u32, opt.spin_mapping.?.quantum_jump_iteration.?, @intCast(j))) {
                    try initialBlochVector(T, &S, &SI, opt.spin_mapping.?, s, rand_bloc);
                }

                if (opt.write.population_mean               != null)  pop.ptr(j, 1 + s).* += 1                                                                  ;
                if (opt.write.kinetic_energy_mean           != null) ekin.ptr(j, 1 + 0).* += Ekin                                                               ;
                if (opt.write.potential_energy_mean         != null) epot.ptr(j, 1 + 0).* += Epot                                                               ;
                if (opt.write.total_energy_mean             != null) etot.ptr(j, 1 + 0).* += Ekin + Epot                                                        ;
                if (opt.write.position_mean                 != null) for (0..r.rows) |k| {position.ptr(j, 1 + k).* += r.at(k);}                                 ;
                if (opt.write.momentum_mean                 != null) for (0..v.rows) |k| {momentum.ptr(j, 1 + k).* += v.at(k) * opt.initial_conditions.mass[k];};
                if (opt.write.time_derivative_coupling_mean != null) for (0..TDC.rows * TDC.cols) |k| {tdc.ptr(j, 1 + k).* += TDC.data[k];}                     ;
                if (opt.write.fssh_coefficient_mean         != null) for (0..C.rows) |k| {coef.ptr(j, 1 + k).* += C.at(k).magnitude() * C.at(k).magnitude();}   ;

                if (j == opt.iterations) assignOutput(T, &output, r, v, s, Ekin, Epot, opt.initial_conditions.mass);

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
        if (print) try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.pop.at(i)});
    }

    if (opt.write.kinetic_energy_mean          ) |path| try     ekin.write(path);
    if (opt.write.population_mean              ) |path| try      pop.write(path);
    if (opt.write.potential_energy_mean        ) |path| try     epot.write(path);
    if (opt.write.total_energy_mean            ) |path| try     etot.write(path);
    if (opt.write.position_mean                ) |path| try position.write(path);
    if (opt.write.momentum_mean                ) |path| try momentum.write(path);
    if (opt.write.time_derivative_coupling_mean) |path| try      tdc.write(path);
    if (opt.write.fssh_coefficient_mean        ) |path| try     coef.write(path);

    return output;
}

/// Diagonalize the potential, assign it to the original potential and correct the sign.
pub fn adiabatizePotential(comptime T: type, U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), UC2: []Matrix(T), T1: *Matrix(T), T2: *Matrix(T), i: usize) void {
    mat.eigh(T, UA, UC, U.*, T1, T2);

    @memcpy(U.data, UA.data); @memcpy(UC2[i % 2].data, UC.data);

    if (i > 0) for (0..UC.cols) |k| {

        var overlap: T = 0;

        for (0..UC.rows) |l| {
            overlap += UC2[i % 2].at(l, k) * UC2[(i - 1) % 2].at(l, k);
        }

        if (overlap < 0) for (0..UC.rows) |l| {
            UC2[i % 2].ptr(l, k).* *= -1;
        };
    };
}

/// Add the position, momentum, population and energies to the output vector.
pub fn assignOutput(comptime T: type, output: *ClassicalDynamicsOutput(T), r: Vector(T), v: Vector(T), state: u32, Ekin: T, Epot: T, mass: []const T) void {
    output.pop.ptr(state).* += 1;

    for (0..r.rows) |i| {
        output.r.ptr(i).* += r.at(i); output.p.ptr(i).* += v.at(i) * mass[i];
    }

    output.Ekin += Ekin;
    output.Epot += Epot;
}

/// Calculate force acting on a specific coordinate "c" multiplied by mass as a negative derivative of the potential.
pub fn calculateForce(comptime T: type, opt: ClassicalDynamicsOptions(T), U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), T1: *Matrix(T), T2: *Matrix(T), r: *Vector(T), c: usize, s: u32) !T {
    r.ptr(c).* += 1 * opt.derivative_step; try mpt.eval(T, U, opt.potential, r.*);

    if (opt.adiabatic) {mat.eigh(T, UA, UC, U.*, T1, T2); @memcpy(U.data, UA.data);} const Up = U.at(s, s);

    r.ptr(c).* -= 2 * opt.derivative_step; try mpt.eval(T, U, opt.potential, r.*);

    if (opt.adiabatic) {mat.eigh(T, UA, UC, U.*, T1, T2); @memcpy(U.data, UA.data);} const Um = U.at(s, s);

    r.ptr(c).* += opt.derivative_step;

    return -0.5 * (Up - Um) / opt.derivative_step;
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
        for (0..FC.rows) |i| {
            K.ptr(i).* = FC.at(i).mul(Complex(T).init(FU.at(i, i), 0)).mulbyi().neg();
        }

        for (0..FC.rows) |i| for (0..FC.rows) |j| {
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
pub fn initialBlochVector(comptime T: type, S: *Matrix(T), SI: *Matrix(usize), opt: ClassicalDynamicsOptions(T).SpinMapping, s: u32, rand: std.Random) !void {
    S.ptr(0, 0).* = 0; S.ptr(0, 1).* = 0; S.ptr(0, 2).* = if (s == 1) 1 else -1;

    for (0..SI.rows) |i| {
        SI.ptr(i, 0).* = s; SI.ptr(i, 1).* = if (i >= s) i + 1 else i;
    }

    // try SI.print(std.io.getStdOut().writer());

    if (!opt.fewest_switches) for (0..S.rows) |i| {

        const phi = 2 * std.math.pi * rand.float(T);

        const cos_theta = std.math.sqrt(rand.float(T));
        const sin_theta = std.math.sqrt(1 - cos_theta * cos_theta);

        S.ptr(i, 0).* = sin_theta * std.math.cos(phi);
        S.ptr(i, 1).* = sin_theta * std.math.sin(phi);
        S.ptr(i, 2).* = cos_theta;
    };
}

/// Function to propagate the vector on the Bloch sphere for the spin mapping methods. The function returns the new state, if a switch occurs.
pub fn spinMapping(comptime T: type, S: *Matrix(T), SI: *Matrix(usize), opt: ClassicalDynamicsOptions(T).SpinMapping, U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, rand: std.Random, SP: *Matrix(T), SN: *Matrix(T)) !u32 {
    var sn = s;

    const OmegaExp = struct { fn get (E: *Matrix(T), VV: T, TT: T) void {
        const a = VV * VV + 4 * TT * TT; const b = Complex(T).init(0, std.math.sqrt(a));
        
        const sinhb = std.math.complex.sinh(b); const coshb = std.math.complex.cosh(b);

        E.ptr(0, 0).* = coshb.re;
        E.ptr(0, 1).* = -VV * sinhb.div(b).re;
        E.ptr(0, 2).* = 2 * TT * sinhb.div(b).re;
        E.ptr(1, 0).* = -E.at(0, 1);
        E.ptr(1, 1).* = (4 * TT * TT + VV * VV * coshb.re) / a;
        E.ptr(1, 2).* = -2 * TT * VV * (coshb.re - 1) / a;
        E.ptr(2, 0).* = -E.at(0, 2);
        E.ptr(2, 1).* = E.at(1, 2);
        E.ptr(2, 2).* = (VV * VV + 4 * TT * TT * coshb.re) / a;
    }};

    for (0..S.rows) |i| {

        const nu = max(SI.at(i, 0), SI.at(i, 1)); const nd = min(SI.at(i, 0), SI.at(i, 1));
        
        OmegaExp.get(SP, (U.at(nu, nu) - U.at(nd, nd)) * time_step, TDC.at(nd, nu) * time_step);

        var SM = S.row(i); var SNM = SN.row(i); std.mem.swap(usize, &SM.rows, &SM.cols); std.mem.swap(usize, &SNM.rows, &SNM.cols);

        mat.mm(T, &SNM, SP.*, SM);
    }

    const rn = rand.float(T); if (opt.fewest_switches) {
        if (s == sn and s == 0 and rn < -2 * TDC.at(0, 1) * SN.at(0, 0) / (1 - SN.at(0, 2)) * time_step) sn = 1;
        if (s == sn and s == 1 and rn <  2 * TDC.at(0, 1) * SN.at(0, 0) / (1 + SN.at(0, 2)) * time_step) sn = 0;
    }

    if (!opt.fewest_switches) {

        for (0..S.rows) |i| if (S.at(i, 2) * SN.at(i, 2) < 0) {
            sn = if (SI.at(i, 0) == s) @intCast(SI.at(i, 1)) else @intCast(SI.at(i, 0));
        };

        for (0..SI.rows) |i| {
            if (SI.at(i, 0) == s and SI.at(i, 1) != sn) SI.ptr(i, 0).* = sn;
            if (SI.at(i, 1) == s and SI.at(i, 0) != sn) SI.ptr(i, 1).* = sn;
        }
    }

    @memcpy(S.data, SN.data); return sn;
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
pub fn printIteration(comptime T: type, i: u32, j: u32, Ekin: T, Epot: T, Etot: T, s: u32, r: Vector(T), v: Vector(T), C: Vector(Complex(T)), S: Matrix(T), mass: []const T, fssh: bool, mash: bool) !void {
    try std.io.getStdOut().writer().print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:5} [", .{i + 1, j, Ekin, Epot, Etot, s});

    for (0..min(r.rows, 3)) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", r.at(k)});
    }

    if (r.rows > 3) {try std.io.getStdOut().writer().print(", ...", .{});}

    try std.io.getStdOut().writer().print("] [", .{});

    for (0..min(v.rows, 3)) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", v.at(k) * mass[k]});
    }

    if (v.rows > 3) {try std.io.getStdOut().writer().print(", ...", .{});}

    if (fssh or mash) try std.io.getStdOut().writer().print("] [", .{});

    if (fssh) for (0..C.rows) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", C.at(k).magnitude() * C.at(k).magnitude()});
    };

    if (mash) for (0..S.cols) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", S.at(0, k)});
    };

    try std.io.getStdOut().writer().print("]\n", .{});
}

/// Function to print the initial header.
pub fn printHeader(ndim: usize, nstate: usize, fssh: bool, mash: bool) !void {
    try std.io.getStdOut().writer().print("\n{s:6} {s:6} {s:12} {s:12} {s:12} {s:5}", .{"TRAJ", "ITER", "EKIN", "EPOT", "ETOT", "STATE"});

    if (ndim > 1) for (0..min(ndim, 3) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});};

    if (ndim > 3) try std.io.getStdOut().writer().print("     ", .{}); try std.io.getStdOut().writer().print(" {s:11}", .{"POSITION"});

    if (ndim > 1) for (0..min(ndim, 3) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});};

    if (ndim > 3) try std.io.getStdOut().writer().print("     ", .{}); try std.io.getStdOut().writer().print(" {s:11}", .{"MOMENTUM"});

    if (fssh) {for (0..min(nstate, 4) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});} try std.io.getStdOut().writer().print(" {s:11}", .{"|COEFS|^2" });}
    if (mash) {for (0..2                 ) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});} try std.io.getStdOut().writer().print(" {s:11}", .{"|BLOCHV|^2"});}

    try std.io.getStdOut().writer().print("\n", .{});
}

/// Function to propagate the classical coordinates in time on an "s" state.
pub fn propagate(comptime T: type, opt: ClassicalDynamicsOptions(T), r: *Vector(T), v: *Vector(T), a: *Vector(T), U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), T1: *Matrix(T), T2: *Matrix(T), s: u32) !void {
    for (0..r.rows) |i| {

        const F = try calculateForce(T, opt, U, UA, UC, T1, T2, r, i, s);

        const ap = a.at(i);

        a.ptr(i).* = F / opt.initial_conditions.mass[i];
        v.ptr(i).* += 0.5 * (a.at(i) + ap) * opt.time_step;
        r.ptr(i).* += (v.at(i) + 0.5 * a.at(i) * opt.time_step) * opt.time_step;
    }
}

/// Function to rescale velocity after a nonadiabatic jump.
pub fn rescaleVelocity(comptime T: type, v: *Vector(T), ps: u32, ns: u32, U: Matrix(T), Ekin: T) void {
    vec.muls(T, v, v.*, std.math.sqrt((Ekin - U.at(ns, ns) + U.at(ps, ps)) / Ekin));
}

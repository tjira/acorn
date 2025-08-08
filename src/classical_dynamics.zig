//! Module for classical dynamics simulations.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const A2AU = @import("constant.zig").A2AU;
const AN2M = @import("constant.zig").AN2M;

const inp = @import("input.zig"    );
const mat = @import("matrix.zig"   );
const pot = @import("potential.zig");
const mth = @import("math.zig"     );
const out = @import("output.zig"   );
const sys = @import("system.zig"   );
const vec = @import("vector.zig"   );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat                  = @import("helper.zig").asfloat                 ;
const extractGeometryFromMovie = @import("helper.zig").extractGeometryFromMovie;

/// Main function for running classical dynamics simulations.
pub fn run(comptime T: type, opt: inp.ClassicalDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !out.ClassicalDynamicsOutput(T) {
    try checkErrors(T, opt, allocator);

    const abinitio = if (opt.hamiltonian.abinitio != null) true else false;

    var potential = try loadPotential(T, opt.hamiltonian, opt.initial_conditions, allocator); defer potential.deinit();

    var mass = try allocator.alloc(T, potential.dims); defer allocator.free(mass); if (!abinitio) @memcpy(mass, opt.initial_conditions.mass.?);

    if (opt.initial_conditions.position_mean != null and opt.initial_conditions.position_mean.?.len != potential.dims) return error.InvalidMeanPosition;
    if (opt.initial_conditions.position_std  != null and opt.initial_conditions.position_std.?.len  != potential.dims) return error.InvalidStdPosition ;
    if (opt.initial_conditions.momentum_mean != null and opt.initial_conditions.momentum_mean.?.len != potential.dims) return error.InvalidMeanMomentum;
    if (opt.initial_conditions.momentum_std  != null and opt.initial_conditions.momentum_std.?.len  != potential.dims) return error.InvalidStdMomentum ;
    if (opt.initial_conditions.mass          != null and opt.initial_conditions.mass.?.len          != potential.dims) return error.InvalidMass        ;

    if (opt.initial_conditions.state.len         != potential.states ) return error.InvalidInitialState;
    if (mth.sum(T, opt.initial_conditions.state) != 1                ) return error.InvalidStateSum    ;
    if (opt.write.bloch_vector != null and potential.states != 2     ) return error.SpinNotCalculable  ;
    if (opt.write.bloch_vector_mean != null and potential.states != 2) return error.SpinNotCalculable  ;

    var output = try out.ClassicalDynamicsOutput(T).init(potential.dims, potential.states, allocator);

    var prng_jump = std.Random.DefaultPrng.init(opt.seed); const rand_jump = prng_jump.random();
    var prng_traj = std.Random.DefaultPrng.init(opt.seed); const rand_traj = prng_traj.random();
    var prng_bloc = std.Random.DefaultPrng.init(opt.seed); const rand_bloc = prng_bloc.random();

    var hopgeoms = std.ArrayList(Vector(T)).init(allocator); defer hopgeoms.deinit();

    var bloch_mean    = try Matrix(T).init(opt.iterations + 1, 1 + 4                                  , allocator); defer    bloch_mean.deinit();
    var coef_mean     = try Matrix(T).init(opt.iterations + 1, 1 + potential.states                   , allocator); defer     coef_mean.deinit();
    var ekin_mean     = try Matrix(T).init(opt.iterations + 1, 1 + 1                                  , allocator); defer     ekin_mean.deinit();
    var epot_mean     = try Matrix(T).init(opt.iterations + 1, 1 + 1                                  , allocator); defer     epot_mean.deinit();
    var etot_mean     = try Matrix(T).init(opt.iterations + 1, 1 + 1                                  , allocator); defer     etot_mean.deinit();
    var momentum_mean = try Matrix(T).init(opt.iterations + 1, 1 + potential.dims                     , allocator); defer momentum_mean.deinit();
    var pop_mean      = try Matrix(T).init(opt.iterations + 1, 1 + potential.states                   , allocator); defer      pop_mean.deinit();
    var position_mean = try Matrix(T).init(opt.iterations + 1, 1 + potential.dims                     , allocator); defer position_mean.deinit();
    var potmat_mean   = try Matrix(T).init(opt.iterations + 1, 1 + potential.states * potential.states, allocator); defer   potmat_mean.deinit();
    var tdc_mean      = try Matrix(T).init(opt.iterations + 1, 1 + potential.states * potential.states, allocator); defer      tdc_mean.deinit();

    bloch_mean   .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    coef_mean    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    ekin_mean    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    epot_mean    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    etot_mean    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    momentum_mean.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    pop_mean     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    position_mean.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    potmat_mean  .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    tdc_mean     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));

    var bloch:    Matrix(T) = undefined;
    var coef:     Matrix(T) = undefined;
    var ekin:     Matrix(T) = undefined;
    var epot:     Matrix(T) = undefined;
    var etot:     Matrix(T) = undefined;
    var momentum: Matrix(T) = undefined;
    var pop:      Matrix(T) = undefined;
    var position: Matrix(T) = undefined;
    var potmat:   Matrix(T) = undefined;
    var tdc:      Matrix(T) = undefined;

    if (opt.write.coefficient != null             ) coef     = try Matrix(T).init(opt.iterations + 1, 1 +                    potential.states * opt.trajectories, allocator);
    if (opt.write.kinetic_energy != null          ) ekin     = try Matrix(T).init(opt.iterations + 1, 1 +                                       opt.trajectories, allocator);
    if (opt.write.potential_energy != null        ) epot     = try Matrix(T).init(opt.iterations + 1, 1 +                                       opt.trajectories, allocator);
    if (opt.write.total_energy != null            ) etot     = try Matrix(T).init(opt.iterations + 1, 1 +                                       opt.trajectories, allocator);
    if (opt.write.momentum != null                ) momentum = try Matrix(T).init(opt.iterations + 1, 1 +                      potential.dims * opt.trajectories, allocator);
    if (opt.write.population != null              ) pop      = try Matrix(T).init(opt.iterations + 1, 1 +                    potential.states * opt.trajectories, allocator);
    if (opt.write.position != null                ) position = try Matrix(T).init(opt.iterations + 1, 1 +                      potential.dims * opt.trajectories, allocator);
    if (opt.write.bloch_vector != null            ) bloch    = try Matrix(T).init(opt.iterations + 1, 1 +                                   4 * opt.trajectories, allocator);
    if (opt.write.time_derivative_coupling != null) tdc      = try Matrix(T).init(opt.iterations + 1, 1 + potential.states * potential.states * opt.trajectories, allocator);
    if (opt.write.potential_matrix != null        ) potmat   = try Matrix(T).init(opt.iterations + 1, 1 + potential.states * potential.states * opt.trajectories, allocator);

    if (opt.write.coefficient              != null) coef    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.kinetic_energy           != null) ekin    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.potential_energy         != null) epot    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.total_energy             != null) etot    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.momentum                 != null) momentum.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.population               != null) pop     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.position                 != null) position.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.bloch_vector             != null) bloch   .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.time_derivative_coupling != null) tdc     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    if (opt.write.potential_matrix         != null) potmat  .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));

    const fssh = opt.fewest_switches != null;
    const lzsh = opt.landau_zener    != null;
    const mash = opt.spin_mapping    != null;

    if ((fssh and lzsh and mash) or (fssh and lzsh) or (fssh and mash) or (lzsh and mash)) {
        return error.MultipleHoppingMechanisms;
    }

    if (mash and potential.states != 2) return error.IncompatibleSurfaceHoppingAlgorithm;

    const tdc_hst    = if (opt.time_derivative_coupling != null) std.mem.eql(u8, opt.time_derivative_coupling.?, "hst"   ) else false;
    const tdc_kappa  = if (opt.time_derivative_coupling != null) std.mem.eql(u8, opt.time_derivative_coupling.?, "kappa" ) else false;
    const tdc_lambda = if (opt.time_derivative_coupling != null) std.mem.eql(u8, opt.time_derivative_coupling.?, "lambda") else false;
    const tdc_npi    = if (opt.time_derivative_coupling != null) std.mem.eql(u8, opt.time_derivative_coupling.?, "npi"   ) else false;
    const tdc_nacv   = if (opt.time_derivative_coupling != null) std.mem.eql(u8, opt.time_derivative_coupling.?, "nacv"  ) else false;

    if (opt.time_derivative_coupling != null and !tdc_hst and !tdc_kappa and !tdc_lambda and !tdc_npi and !tdc_nacv) {
        return error.UnknownTimeDerivativeCoupling;
    }

    {
        var T1 = try Matrix(T).init(potential.states, potential.states, allocator); defer T1.deinit();
        var T2 = try Matrix(T).init(potential.states, potential.states, allocator); defer T2.deinit();

        var r = try Vector(T).init(potential.dims, allocator); defer r.deinit();
        var p = try Vector(T).init(potential.dims, allocator); defer p.deinit();
        var v = try Vector(T).init(potential.dims, allocator); defer v.deinit();
        var a = try Vector(T).init(potential.dims, allocator); defer a.deinit();

        var rp = try Vector(T).init(potential.dims, allocator); defer rp.deinit();
        var pp = try Vector(T).init(potential.dims, allocator); defer pp.deinit();
        var vp = try Vector(T).init(potential.dims, allocator); defer vp.deinit();
        var ap = try Vector(T).init(potential.dims, allocator); defer ap.deinit();

        var U    = try Matrix(T).init(potential.states,                                potential.states, allocator); defer    U.deinit();    U.fill(0);
        var UA   = try Matrix(T).init(potential.states,                                potential.states, allocator); defer   UA.deinit();   UA.fill(0);
        var UC   = try Matrix(T).init(potential.states,                                potential.states, allocator); defer   UC.deinit();   UC.fill(0);
        var UCS  = try Matrix(T).init(potential.states,                                potential.states, allocator); defer  UCS.deinit();  UCS.fill(0);
        var TDC  = try Matrix(T).init(potential.states,                                potential.states, allocator); defer  TDC.deinit();  TDC.fill(0);
        var NACV = try Matrix(T).init(potential.dims * (potential.states * potential.states - potential.states) / 2, 1,      allocator); defer NACV.deinit(); NACV.fill(0);

        var P = try Vector(T                  ).init(potential.states, allocator); defer P.deinit();
        var C = try Vector(std.math.Complex(T)).init(potential.states, allocator); defer C.deinit();
        var S = try Vector(T                  ).init(3,      allocator); defer S.deinit();

        var I   = try Matrix(std.math.Complex(T)).init(potential.states, potential.states, allocator); defer   I.deinit();
        var KC1 = try Vector(std.math.Complex(T)).init(C.rows,                             allocator); defer KC1.deinit();
        var KC2 = try Vector(std.math.Complex(T)).init(C.rows,                             allocator); defer KC2.deinit();
        var KC3 = try Vector(std.math.Complex(T)).init(C.rows,                             allocator); defer KC3.deinit();
        var KC4 = try Vector(std.math.Complex(T)).init(C.rows,                             allocator); defer KC4.deinit();

        var SP = try Matrix(T).init(3, 3, allocator); defer SP.deinit();
        var SN = try Vector(T).init(3,    allocator); defer SN.deinit();
        var S0 = try Vector(T).init(3,    allocator); defer S0.deinit();

        var U3    = [3]Matrix(T){try U.clone(),    try U.clone(), try U.clone()   }; defer    U3[0].deinit(); defer    U3[1].deinit(); defer U3[2].deinit();
        var UC2   = [2]Matrix(T){try U.clone(),    try U.clone()                  }; defer   UC2[0].deinit(); defer   UC2[1].deinit()                      ;
        var NACV2 = [2]Matrix(T){try NACV.clone(), try NACV.clone()               }; defer NACV2[0].deinit(); defer NACV2[1].deinit()                      ;

        if (print) try printHeader(potential.dims, potential.states, fssh, mash);

        for (0..opt.trajectories) |i| {

            if (abinitio) removeAbinitioOutputFiles();

            var s = extractInitialState(T, opt.initial_conditions, opt.trajectories, i);

            if (!abinitio) sampleInitialConditions(T, &r, &p, &v, &a, mass, opt.initial_conditions, rand_traj);

            if (abinitio) try initAbinitio(T, &r, &v, &p, &mass, opt.initial_conditions.position_file, opt.initial_conditions.velocity_file, i, allocator);

            if (fssh or lzsh) {C.fill(std.math.Complex(T).init(0, 0)); C.ptr(s).* = std.math.Complex(T).init(1, 0);}
            if (mash        ) {try initialBlochVector(T, &S, s, rand_bloc);                                        }

            const Wcp: T = if (mash) 2 else 1; const Wpp: T = if (mash) 2 * @abs(S.at(2)) else 1;
            S.memcpy(S0); var S3 = [3]u32{s, s, s}; var ns = s; var Ekin: T = 0; var Epot: T = 0;

            for (0..opt.iterations + 1) |j| {

                const time = opt.time_step * asfloat(T, j); r.memcpy(rp); p.memcpy(pp); v.memcpy(vp); a.memcpy(ap);

                if (j == 0) {
                    try potential.evaluate(&U, r, time); if (!abinitio) try adiabatizePotential(T, &U, &UA, &UC, &UC2, j);
                } else {
                    try propagate(T, potential, opt.time_step, opt.derivative_step, mass, &r, &v, &a, &U, &UA, &UC, &UC2, time, s, j, abinitio, allocator);
                }

                U.memcpy(U3[j % 3]); S3[j % 3] = s;

                Ekin = 0; for (0..v.rows) |k| Ekin += 0.5 * mass[k] * v.at(k) * v.at(k); Epot = U.at(s, s); 

                if (abinitio and tdc_nacv) {
                    try readNacvFromFile(T, &NACV, "NACV.mat", allocator); NACV.memcpy(NACV2[j % 2]);
                }

                if (tdc_hst    and j > 1)     derivativeCouplingHst(   T, &TDC, &UCS, &[_]Matrix(T){UC2[j % 2], UC2[(j - 1) % 2]},                         opt.time_step   );
                if (tdc_kappa  and j > 1)     derivativeCouplingKappa( T, &TDC,       &[_]Matrix(T){U3[j % 3],  U3[(j - 1) % 3],   U3[(j - 2) % 3]},       opt.time_step   );
                if (tdc_lambda and j > 1)     derivativeCouplingLambda(T, &TDC,       &[_]Matrix(T){U3[j % 3],  U3[(j - 1) % 3],   U3[(j - 2) % 3]}, v, a, opt.time_step   );
                if (tdc_npi    and j > 1) try derivativeCouplingNpi(   T, &TDC, &UCS, &UC2, &U, r, potential,                                              opt.time_step, j);
                if (tdc_nacv   and j > 1)     derivativeCouplingNacv(  T, &TDC, &[_]Matrix(T){NACV2[j % 2], NACV2[(j - 1) % 2]}, v, 0                                      );

                if (lzsh and j > 1) ns = try landauZener(T, &C, &P, opt.landau_zener.?, &[_]Matrix(T){U3[j % 3], U3[(j - 1) % 3], U3[(j - 2) % 3]}, s, opt.time_step, rand_jump, &I, &KC1);
                if (fssh and j > 1) ns = try fewestSwitches(T, &C, &P, opt.fewest_switches.?, U, TDC, s, opt.time_step, Ekin, rand_jump, &KC1, &KC2, &KC3, &KC4);
                if (mash and j > 1) ns = try spinMapping(T, &S, U, TDC, s, opt.time_step, Ekin, &SP, &SN);

                if (s != ns and Ekin >= U.at(ns, ns) - U.at(s, s)) {
                    rescaleVelocity(T, &v, s, ns, U, Ekin); s = ns; if (opt.write.hopping_geometries != null) try hopgeoms.append(try r.clone());
                }

                if (mash and opt.spin_mapping.?.resample != null) try resampleBlochVector(T, opt.spin_mapping.?.resample.?, &S, v, vp, s, rand_bloc);

                if (fssh or lzsh) calculateBlocVectorFromCoefficients(T, &S, C);

                if (opt.write.coefficient_mean              != null) for (0..C.rows) |k| {coef_mean.ptr(j, 1 + k).* += C.at(k).magnitude() * C.at(k).magnitude();};
                if (opt.write.kinetic_energy_mean           != null) ekin_mean.ptr(j, 1 + 0).* += Wpp * Ekin                                                      ;
                if (opt.write.momentum_mean                 != null) for (0..v.rows) |k| {momentum_mean.ptr(j, 1 + k).* += Wpp * v.at(k) * mass[k];}              ;
                if (opt.write.population_mean               != null) pop_mean.ptr(j, 1 + s).* += Wpp                                                              ;
                if (opt.write.position_mean                 != null) for (0..r.rows) |k| {position_mean.ptr(j, 1 + k).* += Wpp * r.at(k);}                        ;
                if (opt.write.potential_energy_mean         != null) epot_mean.ptr(j, 1 + 0).* += Wpp * Epot                                                      ;
                if (opt.write.potential_matrix_mean         != null) for (0..U.rows * U.cols) |k| {potmat_mean.ptr(j, 1 + k).* += Wpp * U.data[k];}               ;
                if (opt.write.time_derivative_coupling_mean != null) for (0..TDC.rows * TDC.cols) |k| {tdc_mean.ptr(j, 1 + k).* += Wpp * TDC.data[k];}            ;
                if (opt.write.total_energy_mean             != null) etot_mean.ptr(j, 1 + 0).* += Wpp * (Ekin + Epot)                                             ;

                if (opt.write.coefficient              != null) for (0..C.rows) |k| {coef.ptr(j, 1 + k + i * potential.states).* = C.at(k).magnitude() * C.at(k).magnitude();} ;
                if (opt.write.kinetic_energy           != null) ekin.ptr(j, 1 + i).* = Ekin                                                                                    ;
                if (opt.write.momentum                 != null) for (0..v.rows) |k| {momentum.ptr(j, 1 + k + i * potential.dims).* = v.at(k) * mass[k];}                       ;
                if (opt.write.population               != null) pop.ptr(j, 1 + s + i * potential.states).* = 1                                                                 ;
                if (opt.write.position                 != null) for (0..r.rows) |k| {position.ptr(j, 1 + k + i * potential.dims).* = r.at(k);}                                 ;
                if (opt.write.potential_energy         != null) epot.ptr(j, 1 + i).* = Epot                                                                                    ;
                if (opt.write.potential_matrix         != null) for (0..U.rows * U.cols) |k| {potmat.ptr(j, 1 + k + i * potential.states * potential.states).* = U.data[k];}   ;
                if (opt.write.time_derivative_coupling != null) for (0..TDC.rows * TDC.cols) |k| {tdc.ptr(j, 1 + k + i * potential.states * potential.states).* = TDC.data[k];};
                if (opt.write.total_energy             != null) etot.ptr(j, 1 + i).* = Ekin + Epot                                                                             ;

                if (opt.write.bloch_vector != null) {
                    bloch.ptr(j, 1 + i * 4 + 0).* = S.at(0); bloch.ptr(j, 1 + i * 4 + 1).* = S.at(1); bloch.ptr(j, 1 + i * 4 + 2).* = S.at(2);
                }

                if (opt.write.bloch_vector_mean != null) {
                    bloch_mean.ptr(j, 1 + 0).* += Wcp * S.at(0); bloch_mean.ptr(j, 1 + 1).* += Wcp * S.at(1); bloch_mean.ptr(j, 1 + 2).* += Wpp * (if (mash) mth.sgn(S.at(2)) else S.at(2));
                }

                if (j == opt.iterations) assignOutput(T, &output, r, v, s, Ekin, Epot, mass);

                if (print and (i == 0 or (i + 1) % opt.log_intervals.trajectory == 0) and (j % opt.log_intervals.iteration == 0 or S3[j % 3] != S3[(j - 1) % 3])) {
                    try printIteration(T, @intCast(i), @intCast(j), Ekin, Epot, Ekin + Epot, s, r, v, C, S, mass, fssh, mash);
                }
            }
        }
    }

    if (opt.write.bloch_vector != null) for (0..opt.trajectories) |i| for (0..opt.iterations + 1) |j| {
        bloch.ptr(j, 1 + i * 4 + 3).* = std.math.sqrt(bloch.at(j, 1 + i * 4 + 0) * bloch.at(j, 1 + i * 4 + 0) + bloch.at(j, 1 + i * 4 + 1) * bloch.at(j, 1 + i * 4 + 1));
    };

    if (opt.write.bloch_vector_mean != null) for (0..opt.iterations + 1) |j| {
        bloch_mean.ptr(j, 1 + 3).* = std.math.sqrt(bloch_mean.at(j, 1 + 0) * bloch_mean.at(j, 1 + 0) + bloch_mean.at(j, 1 + 1) * bloch_mean.at(j, 1 + 1));
    };

    for (0..opt.iterations + 1) |i| {for (1..bloch_mean.cols   ) |j|    bloch_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..coef_mean.cols    ) |j|     coef_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..ekin_mean.cols    ) |j|     ekin_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..epot_mean.cols    ) |j|     epot_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..etot_mean.cols    ) |j|     etot_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..momentum_mean.cols) |j| momentum_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..pop_mean.cols     ) |j|      pop_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..position_mean.cols) |j| position_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..potmat_mean.cols  ) |j|   potmat_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}
    for (0..opt.iterations + 1) |i| {for (1..tdc_mean.cols     ) |j|      tdc_mean.ptr(i, j).* /= asfloat(T, opt.trajectories);}

    for (0..output.pop.rows) |i| output.pop.ptr(i).* /= asfloat(T, opt.trajectories);
    for (0..output.r.rows  ) |i|   output.r.ptr(i).* /= asfloat(T, opt.trajectories);
    for (0..output.p.rows  ) |i|   output.p.ptr(i).* /= asfloat(T, opt.trajectories);

    output.Ekin /= asfloat(T, opt.trajectories);
    output.Epot /= asfloat(T, opt.trajectories);

    for (0..potential.states) |i| {
        if (print) try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.pop.at(i)});
    }

    if (opt.write.bloch_vector                 ) |path| try         bloch.write(path);
    if (opt.write.bloch_vector_mean            ) |path| try    bloch_mean.write(path);
    if (opt.write.coefficient                  ) |path| try          coef.write(path);
    if (opt.write.coefficient_mean             ) |path| try     coef_mean.write(path);
    if (opt.write.kinetic_energy               ) |path| try          ekin.write(path);
    if (opt.write.kinetic_energy_mean          ) |path| try     ekin_mean.write(path);
    if (opt.write.momentum                     ) |path| try      momentum.write(path);
    if (opt.write.momentum_mean                ) |path| try momentum_mean.write(path);
    if (opt.write.population                   ) |path| try           pop.write(path);
    if (opt.write.population_mean              ) |path| try      pop_mean.write(path);
    if (opt.write.position                     ) |path| try      position.write(path);
    if (opt.write.position_mean                ) |path| try position_mean.write(path);
    if (opt.write.potential_energy             ) |path| try          epot.write(path);
    if (opt.write.potential_energy_mean        ) |path| try     epot_mean.write(path);
    if (opt.write.potential_matrix             ) |path| try        potmat.write(path);
    if (opt.write.potential_matrix_mean        ) |path| try   potmat_mean.write(path);
    if (opt.write.time_derivative_coupling     ) |path| try           tdc.write(path);
    if (opt.write.time_derivative_coupling_mean) |path| try      tdc_mean.write(path);
    if (opt.write.total_energy                 ) |path| try          etot.write(path);
    if (opt.write.total_energy_mean            ) |path| try     etot_mean.write(path);

    if (opt.write.coefficient              != null)     coef.deinit();
    if (opt.write.kinetic_energy           != null)     ekin.deinit();
    if (opt.write.momentum                 != null) momentum.deinit();
    if (opt.write.population               != null)      pop.deinit();
    if (opt.write.position                 != null) position.deinit();
    if (opt.write.potential_energy         != null)     epot.deinit();
    if (opt.write.potential_matrix         != null)   potmat.deinit();
    if (opt.write.time_derivative_coupling != null)      tdc.deinit();
    if (opt.write.total_energy             != null)     etot.deinit();

    if (opt.write.hopping_geometries) |path| {

        const HG = try Matrix(T).init(hopgeoms.items.len, potential.dims, allocator); defer HG.deinit();

        for (0..hopgeoms.items.len) |i| for (0..potential.dims) |j| {HG.ptr(i, j).* = hopgeoms.items[i].at(j);};

        try HG.write(path);
    }

    for (hopgeoms.items) |*e| e.deinit();

    return output;
}

/// Gets the initial state given the state distribution and the trajectory index.
pub fn extractInitialState(comptime T: type, initial_conditions: inp.ClassicalDynamicsOptions(T).InitialConditions, ntraj: usize, i: usize) u32 {
    var s: u32 = undefined;

    for (0..initial_conditions.state.len) |j| {
        if (asfloat(T, i + 1) / asfloat(T, ntraj) <= mth.sum(T, initial_conditions.state[0..j + 1])) {
            s = @intCast(j); break;
        }
    }

    return s;
}

/// Calculate the bloch vector if given the wavefunction coefficients.
pub fn calculateBlocVectorFromCoefficients(comptime T: type, S: *Vector(T), C: Vector(std.math.Complex(T))) void {
    S.ptr(0).* = 2 * (C.at(0).mul(C.at(1).conjugate())).re;
    S.ptr(1).* = 2 * (C.at(0).mul(C.at(1).conjugate())).im;

    S.ptr(2).* = C.at(1).magnitude() * C.at(1).magnitude() - C.at(0).magnitude() * C.at(0).magnitude();
}

/// Function to sample initial conditions.
pub fn sampleInitialConditions(comptime T: type, r: *Vector(T), p: *Vector(T), v: *Vector(T), a: *Vector(T), mass: []const T, initial_conditions: inp.ClassicalDynamicsOptions(T).InitialConditions, rand: std.Random) void {
    if (initial_conditions.position_mean != null) {for (0..r.rows) |i| r.ptr(i).* = initial_conditions.position_mean.?[i] + initial_conditions.position_std.?[i] * rand.floatNorm(T);}
    if (initial_conditions.momentum_mean != null) {for (0..p.rows) |i| p.ptr(i).* = initial_conditions.momentum_mean.?[i] + initial_conditions.momentum_std.?[i] * rand.floatNorm(T);}

    for (0..v.rows) |i| v.ptr(i).* = p.at(i) / mass[i];

    a.fill(0);
}

/// Function to decides what kind of hamiltonian was provided and returns the potential struct.
pub fn loadPotential(comptime T: type, hamiltonian: inp.ClassicalDynamicsOptions(T).Hamiltonian, initial_conditions: inp.ClassicalDynamicsOptions(T).InitialConditions, allocator: std.mem.Allocator) !pot.Potential(T) {
    var potential_map = try pot.getPotentialMap(T, allocator); defer potential_map.deinit();
    
    if (hamiltonian.name     != null) return                                                                                      potential_map.get(hamiltonian.name.?).?;
    if (hamiltonian.file     != null) return try pot.readPotential(T, hamiltonian.dims.?, hamiltonian.file.?,                                                  allocator);
    if (hamiltonian.matrix   != null) return try pot.getPotential(T, hamiltonian.dims.?, hamiltonian.matrix.?,                                                 allocator);
    if (hamiltonian.abinitio != null) return try pot.initAbinitio(T, initial_conditions.state.len, hamiltonian.abinitio.?, initial_conditions.position_file.?, allocator);

    return error.InvalidHamiltonian;
}

/// Diagonalize the potential, assign it to the original potential and correct the sign.
pub fn adiabatizePotential(comptime T: type, U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), UC2: []Matrix(T), i: usize) !void {
    try cwp.Lapack(T).dsyevd(UA, UC, U.*);

    UA.memcpy(U.*); UC.memcpy(UC2[i % 2]);

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
pub fn assignOutput(comptime T: type, output: *out.ClassicalDynamicsOutput(T), r: Vector(T), v: Vector(T), state: u32, Ekin: T, Epot: T, mass: []const T) void {
    output.pop.ptr(state).* += 1;

    for (0..r.rows) |i| {
        output.r.ptr(i).* += r.at(i); output.p.ptr(i).* += v.at(i) * mass[i];
    }

    output.Ekin += Ekin;
    output.Epot += Epot;
}

/// Calculate force acting on a specific coordinate "c" multiplied by mass as a negative derivative of the potential.
pub fn calculateEnergyBias(comptime T: type, potential: pot.Potential(T), derivative_step: T, U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), r: *Vector(T), t: T, c: usize, s: u32, abinitio: bool, allocator: std.mem.Allocator) !T {
    if (abinitio) {

        const E = try mat.read(T, "ENERGY.mat", allocator); defer E.deinit();
        const G = try mat.read(T, "GRADIENT.mat", allocator); defer G.deinit();

        const F0 = -G.at(s * r.rows + c, 0);
        const F1 = -G.at((s + 1) * r.rows + c, 0);

        const E0 = E.at(s, 0);
        const E1 = E.at(s + 1, 0);

        return 10 * (E1 - E0) * (F1 - F0);
    }

    try potential.evaluate(U, r.*, t); try cwp.Lapack(T).dsyevd(UA, UC, U.*); UA.memcpy(U.*);

    const E0 = U.at(s, s);
    const E1 = U.at(s + 1, s + 1);

    r.ptr(c).* += 1 * derivative_step; try potential.evaluate(U, r.*, t);

    try cwp.Lapack(T).dsyevd(UA, UC, U.*); UA.memcpy(U.*);

    const Up0 = U.at(s, s);
    const Up1 = U.at(s + 1, s + 1);

    r.ptr(c).* -= 2 * derivative_step; try potential.evaluate(U, r.*, t);

    try cwp.Lapack(T).dsyevd(UA, UC, U.*); UA.memcpy(U.*);

    const Um0 = U.at(s, s);
    const Um1 = U.at(s + 1, s + 1);

    r.ptr(c).* += derivative_step;

    const F0 = -0.5 * (Up0 - Um0) / derivative_step;
    const F1 = -0.5 * (Up1 - Um1) / derivative_step;

    const res = 1000 * (E1 - E0) * (F1 - F0);

    return res;
    // return -0.5 * r.at(c);
}

/// Calculate force acting on a specific coordinate "c" multiplied by mass as a negative derivative of the potential.
pub fn calculateForce(comptime T: type, potential: pot.Potential(T), derivative_step: T, U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), r: *Vector(T), t: T, c: usize, s: u32, abinitio: bool, allocator: std.mem.Allocator) !T {
    if (abinitio) {
        const G = try mat.read(T, "GRADIENT.mat", allocator); defer G.deinit(); return -G.at(s * r.rows + c, 0);
    }

    r.ptr(c).* += 1 * derivative_step; try potential.evaluate(U, r.*, t);

    try cwp.Lapack(T).dsyevd(UA, UC, U.*); UA.memcpy(U.*); const Up = U.at(s, s);

    r.ptr(c).* -= 2 * derivative_step; try potential.evaluate(U, r.*, t);

    try cwp.Lapack(T).dsyevd(UA, UC, U.*); UA.memcpy(U.*); const Um = U.at(s, s);

    r.ptr(c).* += derivative_step;

    return -0.5 * (Up - Um) / derivative_step;
}

pub fn checkErrors(comptime T: type, opt: inp.ClassicalDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    if (opt.hamiltonian.name == null and opt.hamiltonian.file == null and opt.hamiltonian.matrix == null and opt.hamiltonian.abinitio == null) return error.InvalidHamiltonian;

    if (opt.hamiltonian.name     != null and (opt.hamiltonian.file != null or opt.hamiltonian.matrix != null or opt.hamiltonian.abinitio != null)) return error.MultipleHamiltonianSources;
    if (opt.hamiltonian.file     != null and (opt.hamiltonian.name != null or opt.hamiltonian.matrix != null or opt.hamiltonian.abinitio != null)) return error.MultipleHamiltonianSources;
    if (opt.hamiltonian.matrix   != null and (opt.hamiltonian.file != null or opt.hamiltonian.name   != null or opt.hamiltonian.abinitio != null)) return error.MultipleHamiltonianSources;
    if (opt.hamiltonian.abinitio != null and (opt.hamiltonian.file != null or opt.hamiltonian.name   != null or opt.hamiltonian.matrix   != null)) return error.MultipleHamiltonianSources;

    if (opt.hamiltonian.matrix   != null and opt.hamiltonian.dims   == null) return error.YouNeedToSpecifyHamiltonianDimensions;
    if (opt.hamiltonian.file     != null and opt.hamiltonian.dims   == null) return error.YouNeedToSpecifyHamiltonianDimensions;
    if (opt.hamiltonian.abinitio != null and opt.hamiltonian.dims   != null) return error.DimsAreAmbiguousForAbinitioHamiltonian;

    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.position_file != null) return error.YouCannotSpecifyInitialPositionFile;
    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.position_mean == null) return error.YouNeedToSpecifyInitialPositionMean;
    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.position_std  == null) return error.YouNeedToSpecifyInitialPositionStd;
    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.velocity_file != null) return error.YouCannotSpecifyInitialVelocityFile;
    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.momentum_mean == null) return error.YouNeedToSpecifyInitialMomentumMean;
    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.momentum_std  == null) return error.YouNeedToSpecifyInitialMomentumStd;
    if (opt.hamiltonian.abinitio == null and opt.initial_conditions.mass          == null) return error.YouNeedToSpecifyCoordinateMasses;

    if (opt.hamiltonian.abinitio != null and opt.initial_conditions.position_file == null) return error.YouHaveToSpecifyInitialPositionFileForAbinitioHamiltonian;
    if (opt.hamiltonian.abinitio != null and opt.initial_conditions.position_mean != null) return error.CannotSpecifyInitialPositionMeanForAbinitioHamiltonian;
    if (opt.hamiltonian.abinitio != null and opt.initial_conditions.position_std  != null) return error.CannotSpecifyInitialPositionStdForAbinitioHamiltonian;
    if (opt.hamiltonian.abinitio != null and opt.initial_conditions.momentum_mean != null) return error.CannotSpecifyInitialMomentumMeanForAbinitioHamiltonian;
    if (opt.hamiltonian.abinitio != null and opt.initial_conditions.momentum_std  != null) return error.CannotSpecifyInitialMomentumStdForAbinitioHamiltonian;
    if (opt.hamiltonian.abinitio != null and opt.initial_conditions.mass          != null) return error.CannotSpecifyCoordinateMassesForAbinitioHamiltonian;

    if (opt.hamiltonian.abinitio == null and opt.time_derivative_coupling != null and  std.mem.eql(u8, opt.time_derivative_coupling.?, "nacv")) return error.IncompatibleTimeDerivativeCoupling;
    if (opt.hamiltonian.abinitio != null and opt.time_derivative_coupling != null and !std.mem.eql(u8, opt.time_derivative_coupling.?, "nacv")) return error.IncompatibleTimeDerivativeCoupling;

    var potential_map = try pot.getPotentialMap(T, allocator); defer potential_map.deinit();

    if (opt.hamiltonian.name != null and !potential_map.contains(opt.hamiltonian.name.?)) return error.InvalidHamiltonianName;
}

/// Calculate the nonadiabatic coupling between two states according to Hammes--Schiffer and Tully.
pub fn derivativeCouplingHst(comptime T: type, TDC: *Matrix(T), UCS: *Matrix(T), UC2: []const Matrix(T), time_step: T) void {
    UCS.fill(0); TDC.fill(0);

    for (0..UCS.rows) |i| for (0..UCS.cols) |j| for (0..UCS.rows) |k| {
        UCS.ptr(i, j).* += UC2[1].at(k, i) * UC2[0].at(k, j);
    };

    for (0..TDC.rows) |i| for (i + 1..TDC.cols) |j| {
        TDC.ptr(i, j).* = (UCS.at(i, j) - UCS.at(j, i)) / (2 * time_step); TDC.ptr(j, i).* = -TDC.at(i, j);
    };
}

/// Calculate the nonadiabatic coupling between two states using the Baeck-An method with dZ/dT=0 approximation.
pub fn derivativeCouplingKappa(comptime T: type, TDC: *Matrix(T), U3: []const Matrix(T), time_step: T) void {
    for (0..TDC.rows) |i| for (i + 1..TDC.cols) |j| {

        const Z0 = @abs(U3[0].at(j, j) - U3[0].at(i, i)); const Z1 = @abs(U3[1].at(j, j) - U3[1].at(i, i)); const Z2 = @abs(U3[2].at(j, j) - U3[2].at(i, i));

        const ddZ0 = (Z0 - 2 * Z1 + Z2) / time_step / time_step;

        const sigma = if (ddZ0 > 1e-12) 0.5 * std.math.sqrt(ddZ0 / Z0) else 0;

        TDC.ptr(i, j).* = sigma; TDC.ptr(j, i).* = -TDC.at(i, j);
    };
}

/// Calculate the nonadiabatic coupling between two states using the Baeck-An method and no approximation in TDC conversion.
pub fn derivativeCouplingLambda(comptime T: type, TDC: *Matrix(T), U3: []const Matrix(T), v: Vector(T), a: Vector(T), time_step: T) void {
    for (0..TDC.rows) |i| for (i + 1..TDC.cols) |j| {

        const Z0 = @abs(U3[0].at(j, j) - U3[0].at(i, i)); const Z1 = @abs(U3[1].at(j, j) - U3[1].at(i, i)); const Z2 = @abs(U3[2].at(j, j) - U3[2].at(i, i));

        const dZ0 = (Z0 - Z1) / time_step; const ddZ0 = (Z0 - 2 * Z1 + Z2) / time_step / time_step;

        const av = a.at(0); const vv = v.at(0);

        const sigma = if (@abs(v.at(0)) > 1e-8 and (ddZ0 - av / vv * dZ0) > 1e-8) 0.5 * std.math.sqrt((ddZ0 - av / vv * dZ0) / Z0) else 0;

        TDC.ptr(i, j).* = sigma; TDC.ptr(j, i).* = -TDC.at(i, j);
    };
}

/// Calculate the time derivative coupling as a dot product between the velocity and the nonadiabatic coupling vector. Only available for ab initio potentials.
pub fn derivativeCouplingNacv(comptime T: type, TDC: *Matrix(T), NACV2: []const Matrix(T), v: Vector(T), step: T) void {
    TDC.fill(0); var l: usize = 0;

    for (0..TDC.rows) |i| for (i + 1..TDC.cols) |j| {

        for (0..v.rows) |k| {

            const NACV = NACV2[1].at(l, 0) + (NACV2[0].at(l, 0) - NACV2[1].at(l, 0)) * step;

            TDC.ptr(i, j).* += v.at(k) * NACV; l += 1;
        }

        TDC.ptr(j, i).* = -TDC.at(i, j);
    };
}

/// Read the calculated NACV from a file.
pub fn readNacvFromFile(comptime T: type, NACV: *Matrix(T), path: []const u8, allocator: std.mem.Allocator) !void {
    const NACV_T = try mat.read(T, path, allocator); defer NACV_T.deinit();

    NACV_T.memcpy(NACV.*);
}

/// Calculate the nonadiabatic coupling between two states using Norm Preserving Interpolation.
pub fn derivativeCouplingNpi(comptime T: type, TDC: *Matrix(T), UCS: *Matrix(T), UC2O: []Matrix(T), U: *Matrix(T), r: Vector(T), potential: pot.Potential(T), time_step: T, iter: usize) !void {
    UCS.fill(0); const UC2 = &[_]Matrix(T){UC2O[iter % 2], UC2O[(iter - 1) % 2]};

    for (0..UCS.rows) |i| for (0..UCS.cols) |j| for (0..UCS.rows) |k| {
        UCS.ptr(i, j).* += UC2[1].at(k, i) * UC2[0].at(k, j);
    };

    for (0..UCS.rows) |l| if (UCS.at(l, l) == 0) {

        try potential.evaluate(U, r, time_step * asfloat(T, iter)); for (U.data) |*e| e.* += 1e-14; try adiabatizePotential(T, U, TDC, UCS, UC2O, iter);

        UCS.fill(0);

        for (0..UCS.rows) |i| for (0..UCS.cols) |j| for (0..UCS.rows) |k| {
            UCS.ptr(i, j).* += UC2[1].at(k, i) * UC2[0].at(k, j);
        };

        break;
    };

    cwp.Eigen(T).logm(TDC, UCS.*); mat.divs(T, TDC, TDC.*, time_step);
}

/// Function to propagate the wavefunction coefficients used in the FSSH method. The function returns the new state, if a switch occurs.
pub fn fewestSwitches(comptime T: type, C: *Vector(std.math.Complex(T)), P: *Vector(T), opt: inp.ClassicalDynamicsOptions(T).FewestSwitches, U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, Ekin: T, rand: std.Random, K1: @TypeOf(C), K2: @TypeOf(C), K3: @TypeOf(C), K4: @TypeOf(C)) !u32 {
    const iters = asfloat(T, opt.quantum_substep); const alpha = if (opt.decoherence_alpha == null) std.math.inf(T) else opt.decoherence_alpha.?; var ns = s; var rn: T = 0;

    const Function = struct { fn get (K: *Vector(std.math.Complex(T)), FC: Vector(std.math.Complex(T)), FU: Matrix(T), FTDC: Matrix(T)) void {
        for (0..FC.rows) |i| {
            K.ptr(i).* = FC.at(i).mul(std.math.Complex(T).init(FU.at(i, i), 0)).mulbyi().neg();
        }

        for (0..FC.rows) |i| for (0..FC.rows) |j| {
            K.ptr(i).* = K.at(i).sub(FC.at(j).mul(std.math.Complex(T).init(FTDC.at(i, j), 0)));
        };
    }};

    for (0..opt.quantum_substep) |_| {

        K1.fill(std.math.Complex(T).init(0, 0)); K2.fill(std.math.Complex(T).init(0, 0)); K3.fill(std.math.Complex(T).init(0, 0)); K4.fill(std.math.Complex(T).init(0, 0)); P.fill(0);

        Function.get(K1, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).add(K1.at(j).mul(std.math.Complex(T).init(time_step / 2 / iters, 0)));
        }

        Function.get(K2, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).sub(K1.at(j).mul(std.math.Complex(T).init(time_step / 2 / iters, 0)));
            C.ptr(j).* = C.at(j).add(K2.at(j).mul(std.math.Complex(T).init(time_step / 2 / iters, 0)));
        }

        Function.get(K3, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).sub(K2.at(j).mul(std.math.Complex(T).init(time_step / 2 / iters, 0)));
            C.ptr(j).* = C.at(j).add(K3.at(j).mul(std.math.Complex(T).init(time_step / 1 / iters, 0)));
        }

        Function.get(K4, C.*, U, TDC);

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).sub(K3.at(j).mul(std.math.Complex(T).init(time_step / iters, 0)));
        }

        for (0..C.rows) |j| {
            C.ptr(j).* = C.at(j).add(K1.at(j).add(K2.at(j).mul(std.math.Complex(T).init(2, 0))).add(K3.at(j).mul(std.math.Complex(T).init(2, 0))).add(K4.at(j)).mul(std.math.Complex(T).init(time_step / 6 / iters, 0)));
        }

        for (0..C.rows) |j| if (j != ns) {
            P.ptr(j).* = 2 * TDC.at(ns, j) * C.at(j).mul(C.at(ns).conjugate()).re / (std.math.pow(T, C.at(ns).magnitude(), 2) + 1e-14) * time_step / iters;
        };

        if (mth.sum(T, P.data) > 1) for (0..P.rows) |i| {P.ptr(i).* /= mth.sum(T, P.data);};

        for (1..P.rows) |i| P.ptr(i).* += P.at(i - 1);

        if (mth.sum(T, P.data) > 0) rn = rand.float(T);

        if (mth.sum(T, P.data) > 0 and s == ns) for (0..P.rows) |i| if (rn > (if (i > 0) P.at(i - 1) else 0) and rn < P.at(i)) {
            ns = @intCast(i);
        };

        for (0..C.rows) |j| if (j != ns) {
            C.ptr(j).* = C.at(j).mul(std.math.Complex(T).init(std.math.exp(-0.5 * @abs(U.at(j, j) - U.at(ns, ns)) * time_step / iters / (1 + alpha / Ekin)), 0));
        };

        var sumc: T = 0; for (0..C.rows) |j| if (j != ns) {sumc += C.at(j).magnitude() * C.at(j).magnitude();};

        if (C.at(ns).magnitude() > 0) {
            C.ptr(ns).* = C.at(ns).mul(std.math.Complex(T).init(std.math.sqrt((1 - sumc) / C.at(ns).magnitude() / C.at(ns).magnitude()), 0));
        }
    }

    return ns;
}

// Transform the coordinate file into the initial position vector.
pub fn initAbinitio(comptime T: type, r: *Vector(T), v: *Vector(T), p: *Vector(T), mass: *[]T, position_file: ?[]const u8, velocity_file: ?[]const u8, i: usize, allocator: std.mem.Allocator) !void {

    if (velocity_file == null) {v.fill(0); p.fill(0);}

    else {

        const velocity_geometry = try std.fmt.allocPrint(allocator, "velocity_{d}.xyz", .{i}); defer allocator.free(velocity_geometry);

        try extractGeometryFromMovie(velocity_geometry, velocity_file.?, i);

        const system = try sys.read(T, velocity_geometry, 0, allocator); defer system.deinit();

        @memcpy(v.data, system.coords.data);

        for (v.data) |*e| e.* /= A2AU;
    }

    const position_geometry = try std.fmt.allocPrint(allocator, "geometry_{d}.xyz", .{i}); defer allocator.free(position_geometry);

    try extractGeometryFromMovie(position_geometry, position_file.?, i);

    const system = try sys.read(T, position_geometry, 0, allocator); defer system.deinit();

    @memcpy(r.data, system.coords.data);

    for (0..system.atoms.rows) |j| for (0..3) |k| {
        mass.*[3 * j + k] = AN2M[@as(usize, @intFromFloat(system.atoms.at(j))) - 1] * 1822.888486;
    };

    for (0..p.rows) |j| {
        p.ptr(j).* = v.at(j) * mass.*[j];
    }
}

/// Function to initialize the initial vector on the Bloch sphere for the spin mapping methods.
pub fn initialBlochVector(comptime T: type, S: *Vector(T), s: u32, rand: std.Random) !void {
    const phi = 2 * std.math.pi * rand.float(T);

    const cos_theta = if (s == 1) rand.float(T) else rand.float(T) - 1;
    const sin_theta = std.math.sqrt(1 - cos_theta * cos_theta);

    S.ptr(0).* = sin_theta * std.math.cos(phi);
    S.ptr(1).* = sin_theta * std.math.sin(phi);
    S.ptr(2).* = cos_theta;
}

/// Removes the files generated by ab initio external programs.
pub fn removeAbinitioOutputFiles() void {
    std.fs.cwd().deleteFile("geometry.xyz.all") catch {};

    std.fs.cwd().deleteFile("ENERGY.mat"  ) catch {};
    std.fs.cwd().deleteFile("GRADIENT.mat") catch {};
    std.fs.cwd().deleteFile("NACV.mat"    ) catch {};

    std.fs.cwd().deleteFile("bagel.out.all") catch {};
    std.fs.cwd().deleteFile("orca.out.all" ) catch {};

    std.fs.cwd().deleteTree(".bagel") catch {};
    std.fs.cwd().deleteTree(".orca" ) catch {};
}

/// Function to resample the Bloch sphere at apecified time points.
pub fn resampleBlochVector(comptime T: type, opt: inp.ClassicalDynamicsOptions(T).SpinMapping.Resample, S: *Vector(T), v: Vector(T), vp: Vector(T), s: u32, rand: std.Random) !void {
    if (@abs(S.at(0)) < opt.threshold and @abs(S.at(1)) < opt.threshold) {
        try initialBlochVector(T, S, s, rand);
    }

    for (0..v.rows) |k| if (opt.reflect and vp.at(k) * v.at(k) < 0) {
        try initialBlochVector(T, S, s, rand); break;
    };
}

/// Function to propagate the vector on the Bloch sphere for the spin mapping methods. The function returns the new state, if a switch occurs.
pub fn spinMapping(comptime T: type, S: *Vector(T), U: Matrix(T), TDC: Matrix(T), s: u32, time_step: T, Ekin: T, SP: *Matrix(T), SN: *Vector(T)) !u32 {
    var sn = s;

    const OmegaExp = struct { fn get (E: *Matrix(T), VV: T, TT: T) void {
        const a = VV * VV + 4 * TT * TT; const b = std.math.Complex(T).init(0, std.math.sqrt(a));
        
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

    // TODO: Check the correctnes of the sign in front of TDC.
    OmegaExp.get(SP, (U.at(1, 1) - U.at(0, 0)) * time_step, -TDC.at(0, 1) * time_step);

    const SM = S.matrix(); var SNM = SN.matrix(); try cwp.Blas(T).dgemm(&SNM, SP.*, false, SM, false);

    if (S.at(2) * SN.at(2) < 0) sn = if (SN.at(2) < 0) 0 else 1;

    if (s != sn and Ekin < U.at(sn, sn) - U.at(s, s)) sn = s;

    SN.memcpy(S.*); return sn;
}

/// Function to calculate the Landau-Zener probability of a transition between two states. The function returns the new state, if a switch occurs.
pub fn landauZener(comptime T: type, C: *Vector(std.math.Complex(T)), P: *Vector(T), opt: inp.ClassicalDynamicsOptions(T).LandauZener, U3: []const Matrix(T), s: u32, time_step: T, rand: std.Random, I: *Matrix(std.math.Complex(T)), T1: *Vector(std.math.Complex(T))) !u32 {
    var ns = s; var rn: T = 0; var maxddZ0: T = 0; P.fill(0);

    for (0..U3[0].rows) |i| if (i != s) {

        const Z0 = @abs(U3[0].at(i, i) - U3[0].at(s, s)); const Z1 = @abs(U3[1].at(i, i) - U3[1].at(s, s)); const Z2 = @abs(U3[2].at(i, i) - U3[2].at(s, s));

        const dZ0 = (Z0 - Z1) / time_step; const dZ1 = (Z1 - Z2) / time_step;

        const ddZ0 = (Z0 - 2 * Z1 + Z2) / time_step / time_step;

        if (dZ0 * dZ1 > 0 or (dZ0 * dZ1 < 0 and ddZ0 < 0) or @abs(ddZ0) < 1e-14) continue;

        const g = Z0; const v = std.math.sqrt(Z0 * ddZ0); const delta: T = 0.25 * std.math.pow(T, g, 2) / v;

        var p = std.math.exp(-2 * std.math.pi * delta); if (std.math.isNan(p)) p = 0;

        if (opt.mode != null) {

            if (std.mem.eql(u8, opt.mode.?, "nearest")) {
                if (@abs(@as(i32, @intCast(i)) - @as(i32, @intCast(s))) == 1) {P.ptr(i).* = p;}
            }

            else if (std.mem.eql(u8, opt.mode.?, "maxprob")) {
                if (mth.sum(T, P.data) < p) {P.fill(0); P.ptr(i).* = p;}
            }

            else if (std.mem.eql(u8, opt.mode.?, "maxdiff")) {
                if (@abs(ddZ0) > maxddZ0) {maxddZ0 = @abs(ddZ0); P.fill(0); P.ptr(i).* = p;}
            }

            else return error.UnknownLandauZenerMode;
        }

        else P.ptr(i).* = p; 

        if (C.rows == 2 and p > 0) {

            I.ptr(0, 0).* = std.math.Complex(T).init(std.math.sqrt(1 - p), 0);
            I.ptr(0, 1).* = std.math.Complex(T).init(0 - std.math.sqrt(p), 0);

            I.ptr(1, 0).* = I.at(0, 1).neg(); I.ptr(1, 1).* = I.at(0, 0);

            const phi: T = -0.25 * std.math.pi + delta * (std.math.log(T, std.math.e, delta) - 1) + std.math.complex.arg(mth.gamma(T, std.math.Complex(T).init(1, -delta)));

            const exp1 = std.math.complex.exp(std.math.Complex(T).init(0, -phi));
            const exp2 = std.math.complex.exp(std.math.Complex(T).init(0,  phi));

            I.ptr(0, 0).* = I.at(0, 0).mul(exp1);
            I.ptr(1, 1).* = I.at(1, 1).mul(exp2);

            var T1M = T1.matrix(); try cwp.Blas(T).zgemm(&T1M, I.*, false, C.matrix(), false); T1.memcpy(C.*);
        }
    };

    if (mth.sum(T, P.data) > 1) for (0..P.rows) |i| {P.ptr(i).* /= mth.sum(T, P.data);};

    for (1..P.rows) |i| P.ptr(i).* += P.at(i - 1);

    if (mth.sum(T, P.data) > 0) rn = rand.float(T);

    if (mth.sum(T, P.data) > 0 and s == ns) for (0..P.rows) |i| if (rn > (if (i > 0) P.at(i - 1) else 0) and rn < P.at(i)) {
        if (C.rows != 2) {std.mem.swap(std.math.Complex(T), C.ptr(i), C.ptr(ns));} ns = @intCast(i); break;
    };

    return ns;
}

/// Function to print the results of a single iteration.
pub fn printIteration(comptime T: type, i: u32, j: u32, Ekin: T, Epot: T, Etot: T, s: u32, r: Vector(T), v: Vector(T), C: Vector(std.math.Complex(T)), S: Vector(T), mass: []const T, fssh: bool, mash: bool) !void {
    try std.io.getStdOut().writer().print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:5} [", .{i + 1, j, Ekin, Epot, Etot, s});

    for (0..mth.min(r.rows, 3)) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", r.at(k)});
    }

    if (r.rows > 3) {try std.io.getStdOut().writer().print(", ...", .{});}

    try std.io.getStdOut().writer().print("] [", .{});

    for (0..mth.min(v.rows, 3)) |k| {
        try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (k == 0) "" else ", ", v.at(k) * mass[k]});
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

/// Function to print the initial header.
pub fn printHeader(ndim: usize, nstate: usize, fssh: bool, mash: bool) !void {
    try std.io.getStdOut().writer().print("\n{s:6} {s:6} {s:12} {s:12} {s:12} {s:5}", .{"TRAJ", "ITER", "EKIN", "EPOT", "ETOT", "STATE"});

    if (ndim > 1) for (0..mth.min(ndim, 3) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});};

    if (ndim > 3) try std.io.getStdOut().writer().print("     ", .{}); try std.io.getStdOut().writer().print(" {s:11}", .{"POSITION"});

    if (ndim > 1) for (0..mth.min(ndim, 3) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});};

    if (ndim > 3) try std.io.getStdOut().writer().print("     ", .{}); try std.io.getStdOut().writer().print(" {s:11}", .{"MOMENTUM"});

    if (fssh) {for (0..mth.min(nstate, 4) - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});} try std.io.getStdOut().writer().print(" {s:11}", .{"|COEFS|^2" });}
    if (mash) {for (0..2                     ) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});} try std.io.getStdOut().writer().print(" {s:11}", .{"|BLOCHV|^2"});}

    try std.io.getStdOut().writer().print("\n", .{});
}

/// Function to propagate the classical coordinates in time on an "s" state.
pub fn propagate(comptime T: type, potential: pot.Potential(T), time_step: T, derivative_step: T, mass: []const T, r: *Vector(T), v: *Vector(T), a: *Vector(T), U: *Matrix(T), UA: *Matrix(T), UC: *Matrix(T), UC2: []Matrix(T), t: T, s: u32, j: usize, abinitio: bool, allocator: std.mem.Allocator) !void {
    for (0..r.rows) |i| {
        r.ptr(i).* += (v.at(i) + 0.5 * a.at(i) * time_step) * time_step;
    }

    try potential.evaluate(U, r.*, t); if (!abinitio) try adiabatizePotential(T, U, UA, UC, UC2, j);

    for (0..r.rows) |i| {

        const F = try calculateForce(T, potential, derivative_step, U, UA, UC, r, t, i, s, abinitio, allocator);

        // F += try calculateEnergyBias(T, potential, derivative_step, U, UA, UC, r, t, i, s, abinitio, allocator);

        const ap = a.at(i);

        a.ptr(i).* = F / mass[i]; v.ptr(i).* += 0.5 * (a.at(i) + ap) * time_step;
    }
}

/// Function to rescale velocity after a nonadiabatic jump.
pub fn rescaleVelocity(comptime T: type, v: *Vector(T), ps: u32, ns: u32, U: Matrix(T), Ekin: T) void {
    vec.muls(T, v, v.*, std.math.sqrt((Ekin - U.at(ns, ns) + U.at(ps, ps)) / Ekin));
}

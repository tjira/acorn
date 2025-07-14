//! Module for quantum dynamics simulations.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const ftr = @import("fouriertransform.zig");
const inp = @import("input.zig"           );
const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );
const out = @import("output.zig"          );
const vec = @import("vector.zig"          );
const wfn = @import("wavefunction.zig"    );

const Matrix       = @import("matrix.zig"      ).Matrix      ;
const Vector       = @import("vector.zig"      ).Vector      ;
const Wavefunction = @import("wavefunction.zig").Wavefunction;

const asfloat = @import("helper.zig").asfloat;

/// The main function to run the quantum dynamics simulation.
pub fn run(comptime T: type, opt: inp.QuantumDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !out.QuantumDynamicsOutput(T) {
    try checkErrors(T, opt, allocator); var potential_map = try mpt.getPotentialMap(T, allocator); defer potential_map.deinit();

    var pot: ?mpt.Potential(T) = null;

    if (opt.hamiltonian.name   != null) pot =                                            potential_map.get(opt.hamiltonian.name.?);
    if (opt.hamiltonian.file   != null) pot = try mpt.readPotential(T, opt.hamiltonian.dims.?, opt.hamiltonian.file.?,  allocator);
    if (opt.hamiltonian.matrix != null) pot = try mpt.getPotential(T, opt.hamiltonian.dims.?, opt.hamiltonian.matrix.?, allocator);

    const ndim = pot.?.dims; const nstate = pot.?.states; const rdim = std.math.pow(usize, opt.grid.points, ndim); defer pot.?.deinit();

    var spectrum_points = std.math.pow(usize, 2, std.math.log2(opt.iterations + 1) + opt.spectrum.nearest_power_of_two); if (opt.spectrum.flip) spectrum_points *= 2;

    var output = try out.QuantumDynamicsOutput(T).init(ndim, nstate, opt.mode[0] + opt.mode[1], allocator);

    var bloch    = try Matrix(T).init(opt.iterations + 1, 1 + 4,      allocator); defer    bloch.deinit(); bloch   .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var pop      = try Matrix(T).init(opt.iterations + 1, 1 + nstate, allocator); defer      pop.deinit(); pop     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var ekin     = try Matrix(T).init(opt.iterations + 1, 1 + 1     , allocator); defer     ekin.deinit(); ekin    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var epot     = try Matrix(T).init(opt.iterations + 1, 1 + 1     , allocator); defer     epot.deinit(); epot    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var etot     = try Matrix(T).init(opt.iterations + 1, 1 + 1     , allocator); defer     etot.deinit(); etot    .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var position = try Matrix(T).init(opt.iterations + 1, 1 + ndim  , allocator); defer position.deinit(); position.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var momentum = try Matrix(T).init(opt.iterations + 1, 1 + ndim  , allocator); defer momentum.deinit(); momentum.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var acf      = try Matrix(T).init(opt.iterations + 1, 1 + 2     , allocator); defer      acf.deinit(); acf     .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
    var acft     = try Matrix(T).init(spectrum_points,    1 + 2     , allocator); defer     acft.deinit();
    var spectrum = try Matrix(T).init(spectrum_points,    1 + 1     , allocator); defer spectrum.deinit();

    if ( opt.spectrum.flip) acft.column(0).linspace(-opt.time_step * asfloat(T, acft.rows / 2), opt.time_step * asfloat(T, acft.rows / 2 - 1));
    if (!opt.spectrum.flip) acft.column(0).linspace(0,                                          opt.time_step * asfloat(T, acft.rows     - 1));

    var wavefunction:       Matrix(T) = undefined; defer if (opt.write.wavefunction       != null                                 )       wavefunction.deinit();
    var density:            Matrix(T) = undefined; defer if (opt.write.density            != null                                 )            density.deinit();
    var bohm_positions:     Matrix(T) = undefined; defer if (opt.write.bohm_position      != null and opt.bohmian_dynamics != null)     bohm_positions.deinit();
    var bohm_position_mean: Matrix(T) = undefined; defer if (opt.write.bohm_position_mean != null and opt.bohmian_dynamics != null) bohm_position_mean.deinit();
    var bohm_momenta:       Matrix(T) = undefined; defer if (opt.write.bohm_momentum      != null and opt.bohmian_dynamics != null)       bohm_momenta.deinit();
    var bohm_momentum_mean: Matrix(T) = undefined; defer if (opt.write.bohm_momentum_mean != null and opt.bohmian_dynamics != null) bohm_momentum_mean.deinit();

    if (opt.write.wavefunction != null                                         ) wavefunction       = try Matrix(T).init(rdim,               ndim + 2 * (opt.iterations + 1) * nstate,       allocator);
    if (opt.write.density      != null                                         ) density            = try Matrix(T).init(rdim,               ndim + 1 * (opt.iterations + 1) * nstate,       allocator);
    if (opt.bohmian_dynamics   != null and opt.write.bohm_position      != null) bohm_positions     = try Matrix(T).init(1 + opt.iterations, 1 + ndim * opt.bohmian_dynamics.?.trajectories, allocator);
    if (opt.bohmian_dynamics   != null and opt.write.bohm_momentum      != null) bohm_momenta       = try Matrix(T).init(1 + opt.iterations, 1 + ndim * opt.bohmian_dynamics.?.trajectories, allocator);
    if (opt.bohmian_dynamics   != null and opt.write.bohm_position_mean != null) bohm_position_mean = try Matrix(T).init(1 + opt.iterations, 1 + ndim,                                       allocator);
    if (opt.bohmian_dynamics   != null and opt.write.bohm_momentum_mean != null) bohm_momentum_mean = try Matrix(T).init(1 + opt.iterations, 1 + ndim,                                       allocator);

    {
        var T1 = try Matrix(std.math.Complex(T)).init(rdim, 1, allocator); defer T1.deinit();

        var rvec = try Matrix(T).init(rdim, ndim, allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(rdim, ndim, allocator); defer kvec.deinit();

        var r = try Vector(T).init(ndim, allocator); defer r.deinit(); r.fill(0);
        var p = try Vector(T).init(ndim, allocator); defer p.deinit(); p.fill(0);

        const dr = std.math.pow(T, (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1), asfloat(T, ndim));

        var W0 = try Wavefunction(T).init(ndim, nstate, opt.grid.points, dr, allocator); defer W0.deinit();
        var W  = try Wavefunction(T).init(ndim, nstate, opt.grid.points, dr, allocator); defer  W.deinit();
        var WA = try Wavefunction(T).init(ndim, nstate, opt.grid.points, dr, allocator); defer WA.deinit();

        var bohm_field:    Vector(T) = undefined; defer if (opt.bohmian_dynamics != null)    bohm_field.deinit();
        var bohm_position: Vector(T) = undefined; defer if (opt.bohmian_dynamics != null) bohm_position.deinit();
        var bohm_momentum: Vector(T) = undefined; defer if (opt.bohmian_dynamics != null) bohm_momentum.deinit();

        if (opt.bohmian_dynamics != null) {

            if (opt.write.bohm_position      != null)     bohm_positions.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
            if (opt.write.bohm_momentum      != null)       bohm_momenta.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
            if (opt.write.bohm_position_mean != null) bohm_position_mean.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
            if (opt.write.bohm_momentum_mean != null) bohm_momentum_mean.column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));

            bohm_field    = try Vector(T).init(rdim,                                       allocator);
            bohm_position = try Vector(T).init(ndim * opt.bohmian_dynamics.?.trajectories, allocator);
            bohm_momentum = try Vector(T).init(ndim * opt.bohmian_dynamics.?.trajectories, allocator);

            for (0..bohm_momentum.rows) |i| bohm_momentum.ptr(i).* = opt.initial_conditions.momentum[i % ndim];
        }

        var WOPT = try allocator.alloc(Wavefunction(T), if (opt.mode[0] > 0) opt.mode[0] - 1 else 0); defer allocator.free(WOPT);

        for (0..WOPT.len) |i| {
            WOPT[i] = try Wavefunction(T).init(ndim, nstate, opt.grid.points, dr, allocator);
        }

        var P = try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator); defer P.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        var VS: [3]std.ArrayList(Matrix(std.math.Complex(T))) = undefined;

        var R = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rdim);
        var K = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rdim);

        for (0..3) |i| {

            VS[i] = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rdim);

            for (0..rdim) |_| {

                try VS[i].append(try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator));

                if (i == 0) {
                    try R.append(try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator));
                    try K.append(try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator));
                }
            }
        }

        var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

        try rgridPotentials(T, &VS, pot.?, opt.hamiltonian.cap, rvec, 0, allocator);

        try rgridPropagators(T, &R, VA, VC,         opt.time_step,                              opt.mode[0] > 0); defer R.deinit();
        kgridPropagators(    T, &K, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, opt.mode[0] > 0); defer K.deinit();

        if (opt.write.wavefunction != null) for (0..rdim) |i| for (0..ndim) |j| {
            wavefunction.ptr(i, j).* = rvec.at(i, j);
        };

        if (opt.write.density != null) for (0..rdim) |i| for (0..ndim) |j| {
            density.ptr(i, j).* = rvec.at(i, j);
        };

        for (0..opt.mode[0] + opt.mode[1]) |i| {

            if (print) try std.io.getStdOut().writer().print("\n{s} TIME DYNAMICS #{d}", .{if (i < opt.mode[0]) "IMAGINARY" else "REAL", if (i < opt.mode[0]) i + 1 else i - opt.mode[0] + 1});

            if (i == opt.mode[0]) {
                try rgridPropagators(T, &R, VA, VC, opt.time_step, false);
            }

            if (i < opt.mode[0] or (opt.mode[0] == 0)) {

                wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.gamma, opt.initial_conditions.state); W.normalize();

                if (opt.initial_conditions.adiabatic) {try wfn.diabatize(T, &WA, W, VC); WA.memcpy(W);}
            }

            if (opt.bohmian_dynamics != null) initBohmianTrajectories(T, &bohm_position, rvec, W, opt.bohmian_dynamics.?.seed);

            W.memcpy(W0);

            if (print) try std.io.getStdOut().writer().print("\n{s:6} {s:12} {s:12} {s:12} {s:13}", .{"ITER", "EKIN", "EPOT", "ETOT",  "ACF"});

            if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"POSITION"  });}
            if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"MOMENTUM"  });}
            if (print) {if (W.nstate > 1) for (0..W.nstate - 1) |_| {try std.io.getStdOut().writer().print(" " ** 10, .{});}; try std.io.getStdOut().writer().print(" {s:10}\n", .{"POPULATION"});}

            for (0..opt.iterations + 1) |j| {

                if (j > 0 and opt.bohmian_dynamics != null) {
                    try propagateBohmianTrajectories(T, &bohm_position, &bohm_momentum, &bohm_field, W, rvec, kvec, opt.time_step, opt.initial_conditions.mass, &T1);
                }

                if (pot.?.tdep and j > 0) {

                    try rgridPotentials(T, &VS, pot.?, opt.hamiltonian.cap, rvec, opt.time_step * asfloat(T, j), allocator);

                    try rgridPropagators(T, &R, VA, VC, opt.time_step, opt.mode[0] > 0);
                }

                if (j > 0) try wfn.propagate(T, &W, R, K, &T1);

                if (i < opt.mode[0]) orthogonalize(T, &W, WOPT, i);

                if (opt.adiabatic) try wfn.adiabatize(T, &WA, W, VC);

                const Ekin = try wfn.ekin(T, W, kvec, opt.initial_conditions.mass, &T1); const Epot: T = wfn.epot(T, W, V);

                wfn.density(T, &P, if (opt.adiabatic) WA else W); wfn.position(T, &r, W, rvec); try wfn.momentum(T, &p, W, kvec, &T1);

                const acfi = W0.overlap(W);

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.population       != null) for (0..nstate) |k| {pop.ptr(j, 1 + k).* = P.at(k, k).re;};
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.position         != null) for (0..ndim) |k| {position.ptr(j, 1 + k).* = r.at(k);}   ;
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.momentum         != null) for (0..ndim) |k| {momentum.ptr(j, 1 + k).* = p.at(k);}   ;
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.kinetic_energy   != null) ekin.ptr(j, 1 + 0).* = Ekin                               ;
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.potential_energy != null) epot.ptr(j, 1 + 0).* = Epot                               ;
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.total_energy     != null) etot.ptr(j, 1 + 0).* = Ekin + Epot                        ;

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.bloch_vector != null) {

                    bloch.ptr(j, 1 + 0).* = 2 * P.at(0, 1).re;
                    bloch.ptr(j, 1 + 1).* = 2 * P.at(0, 1).im;

                    bloch.ptr(j, 1 + 2).* = P.at(1, 1).re - P.at(0, 0).re;

                    bloch.ptr(j, 1 + 3).* = std.math.sqrt(bloch.at(j, 1) * bloch.at(j, 1) + bloch.at(j, 2) * bloch.at(j, 2));
                }

                if (i == opt.mode[0] + opt.mode[1] - 1 and (opt.write.autocorrelation_function != null or opt.write.spectrum != null)) {
                    acf.ptr(j, 1 + 0).* = acfi.re; acf.ptr(j, 1 + 1).* = acfi.im;
                }

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.wavefunction != null) for (0..rdim) |k| for (0..nstate) |l| {
                    wavefunction.ptr(k, ndim + 2 * j * nstate + 2 * l + 0).* = (if (opt.adiabatic) WA else W).data.at(k, l).re;
                    wavefunction.ptr(k, ndim + 2 * j * nstate + 2 * l + 1).* = (if (opt.adiabatic) WA else W).data.at(k, l).im;
                };

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.density != null) for (0..rdim) |k| for (0..nstate) |l| {
                    density.ptr(k, ndim + j * nstate + l).* = (if (opt.adiabatic) WA else W).data.at(k, l).magnitude();
                };

                if (opt.bohmian_dynamics != null) {

                    if (opt.write.bohm_position != null) bohm_position.memcpy(bohm_positions.row(j).vector().slice(1, bohm_positions.cols));
                    if (opt.write.bohm_momentum != null) bohm_momentum.memcpy(  bohm_momenta.row(j).vector().slice(1,   bohm_momenta.cols));

                    if (opt.write.bohm_position_mean != null) for (0..ndim) |k| {

                        var mean_position: T = 0; for (0..opt.bohmian_dynamics.?.trajectories) |l| mean_position += bohm_position.at(l * ndim + k);
                        var mean_momentum: T = 0; for (0..opt.bohmian_dynamics.?.trajectories) |l| mean_momentum += bohm_momentum.at(l * ndim + k);

                        if (opt.write.bohm_position_mean != null) bohm_position_mean.ptr(j, 1 + k).* = mean_position / asfloat(T, opt.bohmian_dynamics.?.trajectories);
                        if (opt.write.bohm_momentum_mean != null) bohm_momentum_mean.ptr(j, 1 + k).* = mean_momentum / asfloat(T, opt.bohmian_dynamics.?.trajectories);
                    };
                }

                if (print and (j % opt.log_intervals.iteration == 0)) try printIteration(T, @intCast(j), Ekin, Epot, r, p, P, acfi);

                if (j == opt.iterations) assignOutput(T, &output, r, p, P, Ekin, Epot, i);

                if (j == opt.iterations and opt.mode[0] > 0 and i < opt.mode[0] - 1) W.memcpy(WOPT[i]);

                if (print and j == opt.iterations and opt.mode[0] + opt.mode[1] == 1) {
                    try std.io.getStdOut().writer().print("\nFINAL ENERGY OF THE PROPAGATED WFN: {d:.6}\n", .{Ekin + Epot});
                }

                if (print and j == opt.iterations and nstate > 1) for (0..nstate) |k| {
                    try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (k == 0) "\n" else "", k, output.P[i].at(k, k).re});
                };
            }
        }

        if (opt.mode[0] + opt.mode[1] > 1 and print) for (0..opt.mode[0] + opt.mode[1]) |i| {
            try std.io.getStdOut().writer().print("{s}FINAL ENERGY OF PROPAGATION #{d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i + 1, output.Ekin[i] + output.Epot[i]});
        };

        for (R.items   ) |*e| e.deinit();
        for (K.items   ) |*e| e.deinit();
        for (V.items   ) |*e| e.deinit();
        for (VA.items  ) |*e| e.deinit();
        for (VC.items  ) |*e| e.deinit();
        for (WOPT      ) |*e| e.deinit();
    }

    if (opt.write.spectrum != null or opt.write.transformed_autocorrelation_function != null) try makeSpectrum(T, opt.spectrum, &acft, &spectrum, acf, allocator);

    if (opt.write.autocorrelation_function            ) |path| try          acf.write(path);
    if (opt.write.bloch_vector                        ) |path| try        bloch.write(path);
    if (opt.write.kinetic_energy                      ) |path| try         ekin.write(path);
    if (opt.write.momentum                            ) |path| try     momentum.write(path);
    if (opt.write.population                          ) |path| try          pop.write(path);
    if (opt.write.position                            ) |path| try     position.write(path);
    if (opt.write.potential_energy                    ) |path| try         epot.write(path);
    if (opt.write.spectrum                            ) |path| try     spectrum.write(path);
    if (opt.write.total_energy                        ) |path| try         etot.write(path);
    if (opt.write.transformed_autocorrelation_function) |path| try         acft.write(path);
    if (opt.write.wavefunction                        ) |path| try wavefunction.write(path);
    if (opt.write.density                             ) |path| try      density.write(path);

    if (opt.bohmian_dynamics != null) {if (opt.write.bohm_momentum     ) |path| try       bohm_momenta.write(path);}
    if (opt.bohmian_dynamics != null) {if (opt.write.bohm_momentum_mean) |path| try bohm_momentum_mean.write(path);}
    if (opt.bohmian_dynamics != null) {if (opt.write.bohm_position     ) |path| try     bohm_positions.write(path);}
    if (opt.bohmian_dynamics != null) {if (opt.write.bohm_position_mean) |path| try bohm_position_mean.write(path);}

    return output;
}

/// Assigns the calculated properties to the output struct.
pub fn assignOutput(comptime T: type, output: *out.QuantumDynamicsOutput(T), r: Vector(T), p: Vector(T), P: Matrix(std.math.Complex(T)), Ekin: T, Epot: T, i: usize) void {
    r.memcpy(output.r[i]);
    p.memcpy(output.p[i]);
    P.memcpy(output.P[i]);

    output.Ekin[i] = Ekin;
    output.Epot[i] = Epot;
}

pub fn checkErrors(comptime T: type, opt: inp.QuantumDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    if (opt.hamiltonian.name == null and  opt.hamiltonian.file == null and opt.hamiltonian.matrix == null) return error.InvalidHamiltonian;

    if (opt.hamiltonian.name   != null and (opt.hamiltonian.file != null or  opt.hamiltonian.matrix != null)) return error.InvalidHamiltonian;
    if (opt.hamiltonian.file   != null and (opt.hamiltonian.name != null or  opt.hamiltonian.matrix != null)) return error.InvalidHamiltonian;
    if (opt.hamiltonian.matrix != null and (opt.hamiltonian.file != null or  opt.hamiltonian.name   != null)) return error.InvalidHamiltonian;

    if (opt.hamiltonian.matrix != null and opt.hamiltonian.dims == null) return error.InvalidHamiltonian;
    if (opt.hamiltonian.file   != null and opt.hamiltonian.dims == null) return error.InvalidHamiltonian;

    var potential_map = try mpt.getPotentialMap(T, allocator); defer potential_map.deinit();

    const nstate = if (opt.hamiltonian.name != null) potential_map.get(opt.hamiltonian.name.?).?.states else opt.hamiltonian.matrix.?.len;
    const ndim   = if (opt.hamiltonian.name != null) potential_map.get(opt.hamiltonian.name.?).?.dims   else opt.hamiltonian.dims.?;

    if (opt.initial_conditions.state > nstate - 1     ) return error.InvalidInitialState     ;
    if (opt.initial_conditions.position.len != ndim   ) return error.InvalidInitialPosition  ;
    if (opt.initial_conditions.momentum.len != ndim   ) return error.InvalidInitialMomentum  ;
    if (opt.write.bloch_vector != null and nstate != 2) return error.CantCalculateBlochVector;
}

/// Initialize bohmian trajectories.
pub fn initBohmianTrajectories(comptime T: type, bohm_position: *Vector(T), rvec: Matrix(T), W: Wavefunction(T), seed: usize) void {
    var prng = std.Random.DefaultPrng.init(seed); const rand = prng.random(); const count = bohm_position.rows / @as(usize, @intCast(W.ndim));

    for (0..count) |i| {

        var cumprob: T = 0; const rn = rand.float(T);

        outer: for (0..rvec.rows) |j| for (0..W.nstate) |k| {

            cumprob += W.at(j, k).magnitude() * W.at(j, k).magnitude() * std.math.pow(T, (rvec.at(1, W.ndim - 1) - rvec.at(0, W.ndim - 1)), asfloat(T, W.ndim));

            if (rn < cumprob) {
                rvec.row(j).vector().memcpy(bohm_position.slice(i * W.ndim, (i + 1) * W.ndim)); break :outer;
            }
        };
    }
}

/// Returns the propagators for the k-space grid.
pub fn kgridPropagators(comptime T: type, K: *std.ArrayList(Matrix(std.math.Complex(T))), nstate: usize, kvec: Matrix(T), time_step: T, mass: T, imaginary: bool) void {
    const unit = std.math.Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    for (0..kvec.rows) |i| for (0..nstate) |j| {

        for(0..kvec.cols) |k| K.items[i].ptr(j, j).* = K.items[i].at(j, j).add(std.math.Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));

        K.items[i].ptr(j, j).* = std.math.complex.exp(K.items[i].at(j, j).mul(std.math.Complex(T).init(-0.5 * time_step / mass, 0)).mul(unit));
    };
}

/// Function to transform the autocorrelation function to a spectrum.
pub fn makeSpectrum(comptime T: type, opt: inp.QuantumDynamicsOptions(T).Spectrum, acft: *Matrix(T), spectrum: *Matrix(T), acf: Matrix(T), allocator: std.mem.Allocator) !void {
    var transform = try Vector(std.math.Complex(T)).init(spectrum.rows,    allocator); defer transform.deinit();
    var frequency = try Matrix(T                  ).init(spectrum.rows, 1, allocator); defer frequency.deinit();

    mpt.kgrid(T, &frequency, 0, asfloat(T, spectrum.rows) * (acf.at(1, 0) - acf.at(0, 0)), @intCast(spectrum.rows));

    for (0..acf.rows) |i| {

        const exp = std.math.exp(-opt.gaussian_window_exponent * acf.at(i, 0) * acf.at(i, 0));

        if (opt.flip) {
            transform.ptr(acft.rows / 2 + i).* = std.math.Complex(T).init(acf.at(i, 1), -acf.at(i, 2)).mul(std.math.Complex(T).init(exp, 0));
            transform.ptr(acft.rows / 2 - i).* = std.math.Complex(T).init(acf.at(i, 1),  acf.at(i, 2)).mul(std.math.Complex(T).init(exp, 0));
        } else {
            transform.ptr(i).* = std.math.Complex(T).init(acf.at(i, 1), -acf.at(i, 2)).mul(std.math.Complex(T).init(exp, 0));
        }
    }

    for (0..acft.rows) |i| {
        acft.ptr(i, 1).* = transform.at(i).re; acft.ptr(i, 2).* = transform.at(i).im;
    }

    try ftr.fftn(f64, transform.data, &[_]usize{@intCast(transform.rows)}, -1);

    for (0..spectrum.rows) |i| {
        spectrum.ptr(i, 0).* = frequency.at(i, 0); spectrum.ptr(i, 1).* = transform.at(i).magnitude();
    }

    for (0..spectrum.rows) |i| for (1..spectrum.rows - i) |j| if (spectrum.at(j - 1, 0) > spectrum.at(j, 0)) for (0..spectrum.cols) |k| {
        std.mem.swap(T, spectrum.ptr(j - 1, k), spectrum.ptr(j, k));
    };
}

/// Propagates the Bohmian trajectories.
pub fn propagateBohmianTrajectories(comptime T: type, bohm_position: *Vector(T), bohm_momentum: *Vector(T), bohm_field: *Vector(T), W: Wavefunction(T), rvec: Matrix(T), kvec: Matrix(T), time_step: T, mass: T, T1: *Matrix(std.math.Complex(T))) !void {
    for (0..W.ndim) |k| {

        bohm_field.fill(0);

        for (0..W.nstate) |i| {

            for (0..W.npoint) |j| T1.ptr(j, 0).* = W.at(j, i);

            try ftr.fftn(f64, T1.data, W.shape, -1);

            for (0..W.npoint) |j| T1.ptr(j, 0).* = T1.at(j, 0).mul(std.math.Complex(T).init(0, kvec.at(j, k)));

            try ftr.fftn(f64, T1.data, W.shape, 1);

            for (0..bohm_field.rows) |j| bohm_field.ptr(j).* += W.at(j, i).conjugate().mul(T1.at(j, 0)).im;
        }

        for (0..bohm_field.rows) |j| {
            var densum: T = 0; for (0..W.nstate) |i| {densum += W.at(j, i).magnitude() * W.at(j, i).magnitude();} bohm_field.ptr(j).* /= (densum + 1e-14);
        }

        for (0..bohm_momentum.rows / W.ndim) |i| {

            var q: T = 0; var f: usize = 0; var c: usize = 0;

            for (0..W.ndim) |j| {

                const qj = std.math.clamp((bohm_position.at(i * W.ndim + j) - rvec.at(0, j)) / W.dr, 0, asfloat(T, W.shape[j] - 1));

                q += @mod(qj, 1) * @mod(qj, 1);

                f += @as(usize, @intFromFloat(@floor(qj))) * std.math.pow(usize, @intCast(W.shape[0]), W.ndim - j - 1);
                c += @as(usize, @intFromFloat(@ceil(qj)))  * std.math.pow(usize, @intCast(W.shape[0]), W.ndim - j - 1);
            }

            bohm_momentum.ptr(i * W.ndim + k).* = bohm_field.at(f) * (1 - std.math.sqrt(q)) + bohm_field.at(c) * std.math.sqrt(q);
        }
    }

    for (0..bohm_position.rows) |i| {
        bohm_position.ptr(i).* += bohm_momentum.at(i) * time_step / mass;
    }
}

/// Returns the potential matrices for each point in the space.
pub fn rgridPotentials(comptime T: type, VS: *[3]std.ArrayList(Matrix(std.math.Complex(T))), pot: mpt.Potential(T), cap: ?[]const u8, rvec: Matrix(T), t: T, allocator: std.mem.Allocator) !void {
    const nstate = pot.states;

    var U   = try Matrix(T).init(nstate, nstate, allocator); defer   U.deinit();
    var UA  = try Matrix(T).init(nstate, nstate, allocator); defer  UA.deinit();
    var UC  = try Matrix(T).init(nstate, nstate, allocator); defer  UC.deinit();
    var UCP = try Matrix(T).init(nstate, nstate, allocator); defer UCP.deinit();

    const capexpr: ?cwp.Expression(T) = if (cap != null) try cwp.Expression(T).init(cap.?, rvec.rows, allocator) else null;

    for (0..rvec.rows) |i| {

        UC.memcpy(UCP); try pot.evaluate(&U, rvec.row(i).vector(), t); try cwp.Lapack(T).dsyevd(&UA, &UC, U);

        if (i > 0) for (0..UC.cols) |j| {

            var overlap: T = 0;

            for (0..UC.rows) |k| overlap += UC.at(k, j) * UCP.at(k, j);

            if (overlap < 0) for (0..UC.rows) |k| {UC.ptr(k, j).* *= -1;};
        };

        var UAC = try UA.complex();

        if (capexpr != null) {
            const capv = capexpr.?.evaluate(rvec.row(i).vector(), t); for (0..U.rows) |j| {
                UAC.ptr(j, j).* = UAC.at(j, j).add(std.math.Complex(T).init(0, capv));
            }
        }

        UAC.memcpy(VS[1].items[i]);

        for (0..nstate) |j| for (0..nstate) |k| {
            VS[0].items[i].ptr(j, k).* = std.math.Complex(T).init(U .at(j, k), 0);
            VS[2].items[i].ptr(j, k).* = std.math.Complex(T).init(UC.at(j, k), 0);
        };
    }

    if (capexpr != null) capexpr.?.deinit();
}

/// Returns the propagators for the r-space grid.
pub fn rgridPropagators(comptime T: type, R: *std.ArrayList(Matrix(std.math.Complex(T))), VA: std.meta.Child(@TypeOf(R)), VC: @TypeOf(VA), time_step: T, imaginary: bool) !void {
    const unit = std.math.Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1); var T1 = try R.items[0].clone();

    for (0..R.items.len) |i| {

        R.items[i].fill(std.math.Complex(T).init(0, 0)); for (0..R.items[i].rows) |j| {
            R.items[i].ptr(j, j).* = std.math.complex.exp(VA.items[i].at(j, j).mul(std.math.Complex(T).init(-0.5 * time_step, 0)).mul(unit));
        }

        try cwp.Blas(T).zgemm(&T1, VC.items[i], false, R.items[i], false); try cwp.Blas(T).zgemm(&R.items[i], T1, false, VC.items[i], true);
    }
}

/// Prints the iteration information.
pub fn printIteration(comptime T: type, i: usize, Ekin: T, Epot: T, r: Vector(T), p: Vector(T), P: Matrix(std.math.Complex(T)), acfi: std.math.Complex(T)) !void {
        try std.io.getStdOut().writer().print("{d:6} {d:12.6} {d:12.6} {d:12.6} {d:6.3}{s}{d:5.3}i [", .{i, Ekin, Epot, Ekin + Epot, acfi.re, if (acfi.im < 0) "-" else "+", @abs(acfi.im)});

        for (0..r.rows) |j| {
            try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (j == 0) "" else ", ", r.at(j)});
        }

        try std.io.getStdOut().writer().print("] [", .{});

        for (0..p.rows) |j| {
            try std.io.getStdOut().writer().print("{s}{d:9.4}", .{if (j == 0) "" else ", ", p.at(j)});
        }

        try std.io.getStdOut().writer().print("] [", .{});

        for (0..P.rows) |j| {
            try std.io.getStdOut().writer().print("{s}{d:8.5}", .{if (j == 0) "" else ", ", P.at(j, j).re});
        }

        try std.io.getStdOut().writer().print("]\n", .{});
}

/// Orthogonalizes the wavefunction to the lower eigenstates.
pub fn orthogonalize(comptime T: type, W: *Wavefunction(T), WOPT: []Wavefunction(T), states: usize) void {
    for (0..states) |k| {

        const overlap = WOPT[k].overlap(W.*);

        for (0..W.npoint) |l| for (0..W.nstate) |m| {
            W.ptr(l, m).* = W.at(l, m).sub(WOPT[k].data.at(l, m).mul(overlap));
        };
    }

    W.normalize();
}

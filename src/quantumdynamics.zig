//! Module for quantum dynamics simulations.

const std = @import("std");

const bls = @import("blas.zig"          );
const ftw = @import("fftw.zig"          );
const inp = @import("input.zig"         );
const lpk = @import("lapack.zig"        );
const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");
const out = @import("output.zig"        );
const vec = @import("vector.zig"        );
const wfn = @import("wavefunction.zig"  );

const Matrix       = @import("matrix.zig"      ).Matrix      ;
const Vector       = @import("vector.zig"      ).Vector      ;
const Wavefunction = @import("wavefunction.zig").Wavefunction;

const asfloat = @import("helper.zig").asfloat;

/// The main function to run the quantum dynamics simulation.
pub fn run(comptime T: type, opt: inp.QuantumDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !out.QuantumDynamicsOutput(T) {
    if (opt.hamiltonian.name == null and (opt.hamiltonian.dims == null and opt.hamiltonian.matrix == null)) return error.InvalidHamiltonian;
    if (opt.hamiltonian.name != null and (opt.hamiltonian.dims != null or  opt.hamiltonian.matrix != null)) return error.InvalidHamiltonian;

    var pot = if (opt.hamiltonian.name != null) (try mpt.getPotentialMap(T, allocator)).get(opt.hamiltonian.name.?) else try mpt.getPotential(T, opt.hamiltonian.dims.?, opt.hamiltonian.matrix.?, allocator);

    const ndim = pot.?.dims; const nstate = pot.?.states; const rdim = std.math.pow(u32, opt.grid.points, ndim); defer pot.?.deinit();

    if (pot == null                                   ) return error.UnknownPotential        ;
    if (opt.initial_conditions.state > nstate - 1     ) return error.InvalidInitialState     ;
    if (opt.initial_conditions.position.len != ndim   ) return error.InvalidInitialPosition  ;
    if (opt.initial_conditions.momentum.len != ndim   ) return error.InvalidInitialMomentum  ;
    if (opt.write.bloch_vector != null and nstate != 2) return error.CantCalculateBlochVector;

    var spectrum_points = std.math.pow(u32, 2, std.math.log2(opt.iterations + 1) + opt.spectrum.nearest_power_of_two); if (opt.spectrum.flip) spectrum_points *= 2;

    var output = try out.QuantumDynamicsOutput(T).init(ndim, nstate, opt.mode[0] + opt.mode[1], allocator);

    var bloch    = try Matrix(T).init(opt.iterations + 1, 1 + 3,      allocator); defer    bloch.deinit(); bloch   .column(0).linspace(0, opt.time_step * asfloat(T, opt.iterations));
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
    var bohm_positions:     Matrix(T) = undefined; defer if (opt.write.bohm_position      != null and opt.bohmian_dynamics != null)     bohm_positions.deinit();
    var bohm_position_mean: Matrix(T) = undefined; defer if (opt.write.bohm_position_mean != null and opt.bohmian_dynamics != null) bohm_position_mean.deinit();
    var bohm_momenta:       Matrix(T) = undefined; defer if (opt.write.bohm_momentum      != null and opt.bohmian_dynamics != null)       bohm_momenta.deinit();
    var bohm_momentum_mean: Matrix(T) = undefined; defer if (opt.write.bohm_momentum_mean != null and opt.bohmian_dynamics != null) bohm_momentum_mean.deinit();

    if (opt.write.wavefunction != null                                         ) wavefunction       = try Matrix(T).init(rdim,               ndim + 2 * (opt.iterations + 1) * nstate,       allocator);
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

        const VS = try rgridPotentials(T, pot.?, rvec, allocator); var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

        var R = try rgridPropagators(T, VA, VC,   rvec, opt.time_step,                              opt.mode[0] > 0, allocator); defer R.deinit();
        var K = try kgridPropagators(T, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, opt.mode[0] > 0, allocator); defer K.deinit();

        if (opt.write.wavefunction != null) for (0..rdim) |i| for (0..ndim) |j| {
            wavefunction.ptr(i, j).* = rvec.at(i, j);
        };

        for (0..opt.mode[0] + opt.mode[1]) |i| {

            if (print) try std.io.getStdOut().writer().print("\n{s} TIME DYNAMICS #{d}", .{if (i < opt.mode[0]) "IMAGINARY" else "REAL", if (i < opt.mode[0]) i + 1 else i - opt.mode[0] + 1});

            if (i == opt.mode[0]) {

                for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} R.deinit(); K.deinit();

                R = try rgridPropagators(T, VA, VC,   rvec, opt.time_step,                              false, allocator);
                K = try kgridPropagators(T, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, false, allocator);
            }

            if (i < opt.mode[0] or (opt.mode[0] == 0)) {
                wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.gamma, opt.initial_conditions.state); W.normalize();
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

                if (j > 0) try wfn.propagate(T, &W, R, K, &T1);

                if (i < opt.mode[0]) orthogonalize(T, &W, WOPT, i);

                if (opt.adiabatic) wfn.adiabatize(T, &WA, W, VC);

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
                }

                if (i == opt.mode[0] + opt.mode[1] - 1 and (opt.write.autocorrelation_function != null or opt.write.spectrum != null)) {
                    acf.ptr(j, 1 + 0).* = acfi.re; acf.ptr(j, 1 + 1).* = acfi.im;
                }

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.wavefunction != null) for (0..rdim) |k| for (0..nstate) |l| {
                    wavefunction.ptr(k, ndim + 2 * j * nstate + 2 * l + 0).* = (if (opt.adiabatic) WA else W).data.at(k, l).re;
                    wavefunction.ptr(k, ndim + 2 * j * nstate + 2 * l + 1).* = (if (opt.adiabatic) WA else W).data.at(k, l).im;
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
pub fn kgridPropagators(comptime T: type, nstate: usize, kvec: Matrix(T), time_step: T, mass: T, imaginary: bool, allocator: std.mem.Allocator) !std.ArrayList(Matrix(std.math.Complex(T))) {
    const unit = std.math.Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var K = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, kvec.rows);

    for (0..kvec.rows) |i| {

        try K.append(try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator));

        for (0..nstate) |j| {

            for(0..kvec.cols) |k| K.items[i].ptr(j, j).* = K.items[i].at(j, j).add(std.math.Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));

            K.items[i].ptr(j, j).* = std.math.complex.exp(K.items[i].at(j, j).mul(std.math.Complex(T).init(-0.5 * time_step / mass, 0)).mul(unit));
        }
    }

    return K;
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

    ftw.fftwnd(transform.data, &[_]u32{@intCast(transform.rows)}, -1);

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

            ftw.fftwnd(T1.data, W.shape, -1);

            for (0..W.npoint) |j| T1.ptr(j, 0).* = T1.at(j, 0).mul(std.math.Complex(T).init(0, kvec.at(j, k)));

            ftw.fftwnd(T1.data, W.shape, 1);

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
pub fn rgridPotentials(comptime T: type, pot: mpt.Potential(T), rvec: Matrix(T), allocator: std.mem.Allocator) ![3]std.ArrayList(Matrix(std.math.Complex(T))) {
    const nstate = pot.states;

    var U   = try Matrix(T).init(nstate, nstate, allocator); defer   U.deinit();
    var UA  = try Matrix(T).init(nstate, nstate, allocator); defer  UA.deinit();
    var UC  = try Matrix(T).init(nstate, nstate, allocator); defer  UC.deinit();
    var UCP = try Matrix(T).init(nstate, nstate, allocator); defer UCP.deinit();

    var V  = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rvec.rows);
    var VA = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rvec.rows);
    var VC = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        UC.memcpy(UCP); pot.evaluate(&U, rvec.row(i).vector()); lpk.dsyevd(&UA, &UC, U);

        if (i > 0) for (0..UC.cols) |j| {

            var overlap: T = 0;

            for (0..UC.rows) |k| overlap += UC.at(k, j) * UCP.at(k, j);

            if (overlap < 0) for (0..UC.rows) |k| {UC.ptr(k, j).* *= -1;};
        };

        try V.append(try U.complex()); try VA.append(try UA.complex()); try VC.append(try UC.complex());
    }

    return .{V, VA, VC};
}

/// Returns the propagators for the r-space grid.
pub fn rgridPropagators(comptime T: type, VA: std.ArrayList(Matrix(std.math.Complex(T))), VC: @TypeOf(VA), rvec: Matrix(T), time_step: T, imaginary: bool, allocator: std.mem.Allocator) !@TypeOf(VA) {
    const unit = std.math.Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var T1 = try Matrix(std.math.Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer T1.deinit();
    var T2 = try Matrix(std.math.Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer T2.deinit();

    var R = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        T1.fill(std.math.Complex(T).init(0, 0)); for (0..T1.rows) |j| {
            T1.ptr(j, j).* = std.math.complex.exp(VA.items[i].at(j, j).mul(std.math.Complex(T).init(-0.5 * time_step, 0)).mul(unit));
        }

        bls.zgemm(&T2, VC.items[i], false, T1, false); bls.zgemm(&T1, T2, false, VC.items[i], true); try R.append(try T1.clone());
    }

    return R;
}

/// Prints the iteration information.
pub fn printIteration(comptime T: type, i: u32, Ekin: T, Epot: T, r: Vector(T), p: Vector(T), P: Matrix(std.math.Complex(T)), acfi: std.math.Complex(T)) !void {
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

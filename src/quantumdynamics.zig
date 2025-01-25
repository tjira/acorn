//! Module for quantum dynamics simulations.

const std = @import("std"); const Complex = std.math.Complex;

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );
const vec = @import("vector.zig"          );
const wfn = @import("wavefunction.zig"    );

const Matrix       = @import("matrix.zig"      ).Matrix      ;
const Vector       = @import("vector.zig"      ).Vector      ;
const Wavefunction = @import("wavefunction.zig").Wavefunction;

const asfloat = @import("helper.zig").asfloat;

/// The quantum dynamics options struct.
pub fn QuantumDynamicsOptions(comptime T: type) type {
    return struct {
        pub const Grid = struct {
            limits: []const T, points: u32
        };
        pub const InitialConditions = struct {
            position: []const T, momentum: []const T, gamma: T, state: u32, mass: T
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Write = struct {
            autocorrelation_function: ?[]const u8 = null,
            kinetic_energy: ?[]const u8 = null,
            momentum: ?[]const u8 = null,
            population: ?[]const u8 = null,
            position: ?[]const u8 = null,
            potential_energy: ?[]const u8 = null,
            spectrum: ?[]const u8 = null,
            total_energy: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null
        };

        adiabatic: bool,
        iterations: u32,
        mode: []const u32,
        potential: []const u8,
        time_step: T,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals = .{}, write: Write = .{}
    };
}

/// The quantum dynamics output struct.
pub fn QuantumDynamicsOutput(comptime T: type) type {
    return struct {
        P: []Matrix(T),
        r: []Vector(T),
        p: []Vector(T),

        Ekin: []T,
        Epot: []T,

        allocator: std.mem.Allocator,

        /// Initialize the quantum dynamics output struct.
        pub fn init(ndim: usize, nstate: usize, propagations: usize, allocator: std.mem.Allocator) !QuantumDynamicsOutput(T) {
            var output = QuantumDynamicsOutput(T){
                .P = try allocator.alloc(Matrix(T), propagations),
                .r = try allocator.alloc(Vector(T), propagations),
                .p = try allocator.alloc(Vector(T), propagations),

                .Ekin = try allocator.alloc(T, propagations),
                .Epot = try allocator.alloc(T, propagations),

                .allocator = allocator,
            };

            for (0..propagations) |i| {
                output.P[i] = try Matrix(T).init(nstate, nstate, allocator);
                output.r[i] = try Vector(T).init(ndim,           allocator);
                output.p[i] = try Vector(T).init(ndim,           allocator);
            }

            return output;
        }

        /// Free the memory allocated for the quantum dynamics output struct.
        pub fn deinit(self: QuantumDynamicsOutput(T)) void {
            self.allocator.free(self.Ekin); self.allocator.free(self.Epot);

            for (0..self.P.len) |i| {
                self.P[i].deinit(); self.r[i].deinit(); self.p[i].deinit();
            }

            self.allocator.free(self.P); self.allocator.free(self.r); self.allocator.free(self.p);
        }
    };
}

/// The main function to run the quantum dynamics simulation.
pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !QuantumDynamicsOutput(T) {
    if (opt.initial_conditions.state > try mpt.states(opt.potential) - 1) return error.InvalidInitialState;

    const ndim = try mpt.dims(opt.potential); const nstate = try mpt.states(opt.potential); const rdim = std.math.pow(u32, opt.grid.points, ndim);

    var output = try QuantumDynamicsOutput(T).init(ndim, nstate, opt.mode[0] + opt.mode[1], allocator);

    var pop      = try Matrix(T).init(opt.iterations + 1,                                      1 + nstate, allocator); defer      pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations + 1,                                      1 + 1     , allocator); defer     ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations + 1,                                      1 + 1     , allocator); defer     epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations + 1,                                      1 + 1     , allocator); defer     etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations + 1,                                      1 + ndim  , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations + 1,                                      1 + ndim  , allocator); defer momentum.deinit(); momentum.fill(0);
    var acf      = try Matrix(T).init(opt.iterations + 1,                                      1 + 2     , allocator); defer      acf.deinit();      acf.fill(0);
    var spectrum = try Matrix(T).init(2 * std.math.pow(usize, 2, std.math.log2(acf.rows) + 1), 1 + 1     , allocator); defer spectrum.deinit(); spectrum.fill(0);

    {
        const time = try Matrix(T).init(opt.iterations + 1, 1, allocator); defer time.deinit(); time.linspace(0, opt.time_step * asfloat(T, opt.iterations));

        for (0..opt.iterations + 1) |i| {
            pop.ptr(i, 0).*      = time.at(i, 0);
            ekin.ptr(i, 0).*     = time.at(i, 0);
            epot.ptr(i, 0).*     = time.at(i, 0);
            etot.ptr(i, 0).*     = time.at(i, 0);
            position.ptr(i, 0).* = time.at(i, 0);
            momentum.ptr(i, 0).* = time.at(i, 0);
            acf.ptr(i, 0).*      = time.at(i, 0);
            spectrum.ptr(i, 0).* = time.at(i, 0);
        }
    }

    const wrows = if (opt.write.wavefunction != null) rdim else 0; const wcols = if (opt.write.wavefunction != null) ndim + 2 * (opt.iterations + 1) * nstate else 0;

    var wavefunction = try Matrix(T).init(wrows, wcols, allocator); defer wavefunction.deinit(); wavefunction.fill(0);

    {
        var T1 = try Matrix(Complex(T)).init(rdim, 1, allocator); defer T1.deinit();

        var rvec = try Matrix(T).init(rdim, ndim, allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(rdim, ndim, allocator); defer kvec.deinit();

        var r = try Vector(T).init(ndim, allocator); defer r.deinit(); r.fill(0);
        var p = try Vector(T).init(ndim, allocator); defer p.deinit(); p.fill(0);

        const dr = (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1);

        var W0 = try Wavefunction(T).init(ndim, nstate, opt.grid.points, allocator); defer W0.deinit();
        var W  = try Wavefunction(T).init(ndim, nstate, opt.grid.points, allocator); defer  W.deinit();
        var WA = try Wavefunction(T).init(ndim, nstate, opt.grid.points, allocator); defer WA.deinit();

        var WOPT = try std.ArrayList(Wavefunction(T)).initCapacity(allocator, if (opt.mode[0] > 0) opt.mode[0] - 1 else 0); defer WOPT.deinit();

        for (0..WOPT.capacity) |_| {
            try WOPT.append(try Wavefunction(T).init(ndim, nstate, opt.grid.points, allocator));
        }

        var P = try Matrix(T).init(nstate, nstate, allocator); defer P.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        const VS = try rgridPotentials(T, opt.potential, rvec, allocator); var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

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
                wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.gamma, opt.initial_conditions.state); wfn.normalize(T, &W, dr);
            }

            @memcpy(W0.data.data, W.data.data);

            if (print) try std.io.getStdOut().writer().print("\n{s:6} {s:12} {s:12} {s:12} {s:13}", .{"ITER", "EKIN", "EPOT", "ETOT",  "ACF"});

            if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"POSITION"  });}
            if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"MOMENTUM"  });}
            if (print) {if (W.nstate > 1) for (0..W.nstate - 1) |_| {try std.io.getStdOut().writer().print(" " ** 10, .{});}; try std.io.getStdOut().writer().print(" {s:10}\n", .{"POPULATION"});}

            for (0..opt.iterations + 1) |j| {

                if (j > 0) try wfn.propagate(T, &W, R, K, &T1);

                if (i < opt.mode[0]) {
                    
                    for (0..i) |k| {

                        const overlap = wfn.overlap(T, WOPT.items[k], W, dr);

                        for (0..rdim) |l| for (0..nstate) |m| {
                            W.data.ptr(l, m).* = W.data.at(l, m).sub(WOPT.items[k].data.at(l, m).conjugate().mul(overlap));
                        };
                    }

                    wfn.normalize(T, &W, dr);
                }

                if (opt.adiabatic) wfn.adiabatize(T, &WA, W, VC);

                const Ekin = try wfn.ekin(T, W, kvec, opt.initial_conditions.mass, dr, &T1); const Epot: T = wfn.epot(T, W, V, dr);

                wfn.density(T, &P, if (opt.adiabatic) WA else W, dr); wfn.position(T, &r, W, rvec, dr); try wfn.momentum(T, &p, W, kvec, dr, &T1);

                const acfi = wfn.overlap(T, W0, W, dr);

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.population       != null) for (0..nstate) |k| {pop.ptr(j, 1 + k).* = P.at(k, k);};
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.position         != null) for (0..ndim) |k| {position.ptr(j, 1 + k).* = r.at(k);};
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.momentum         != null) for (0..ndim) |k| {momentum.ptr(j, 1 + k).* = p.at(k);};
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.kinetic_energy   != null) ekin.ptr(j, 1 + 0).* = Ekin                            ;
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.potential_energy != null) epot.ptr(j, 1 + 0).* = Epot                            ;
                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.total_energy     != null) etot.ptr(j, 1 + 0).* = Ekin + Epot                     ;

                if (i == opt.mode[0] + opt.mode[1] - 1 and (opt.write.autocorrelation_function != null or opt.write.spectrum != null)) {
                    acf.ptr(j, 1 + 0).* = acfi.re; acf.ptr(j, 1 + 1).* = acfi.im;
                }

                if (i == opt.mode[0] + opt.mode[1] - 1 and opt.write.wavefunction != null) for (0..rdim) |k| for (0..nstate) |l| {
                    wavefunction.ptr(k, ndim + 2 * j * nstate + 2 * l + 0).* = (if (opt.adiabatic) WA else W).data.at(k, l).re;
                    wavefunction.ptr(k, ndim + 2 * j * nstate + 2 * l + 1).* = (if (opt.adiabatic) WA else W).data.at(k, l).im;
                };

                if (print and (j % opt.log_intervals.iteration == 0)) try printIteration(T, @intCast(j), Ekin, Epot, r, p, P, acfi);

                if (j == opt.iterations) {

                    if (opt.mode[0] > 0 and i < opt.mode[0] - 1) @memcpy(WOPT.items[i].data.data, W.data.data);

                    @memcpy(output.P[i].data, P.data); @memcpy(output.r[i].data, r.data); @memcpy(output.p[i].data, p.data);

                    output.Ekin[i] = Ekin; output.Epot[i] = Epot;

                    if (opt.mode[0] + opt.mode[1] == 1 and print) try std.io.getStdOut().writer().print("\nFINAL ENERGY OF THE PROPAGATED WFN: {d:.6}\n", .{Ekin + Epot});

                    if (nstate > 1) for (0..nstate) |k| if (print) {
                        try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (k == 0) "\n" else "", k, output.P[i].at(k, k)});
                    };
                }
            }
        }

        if (opt.mode[0] + opt.mode[1] > 1 and print) for (0..opt.mode[0] + opt.mode[1]) |i| {
            try std.io.getStdOut().writer().print("{s}FINAL ENERGY OF PROPAGATION #{d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i + 1, output.Ekin[i] + output.Epot[i]});
        };

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();} for (VA.items) |*e| {e.deinit();} for (VC.items) |*e| {e.deinit();} for (WOPT.items) |*e| {e.deinit();}
    }

    if (opt.write.spectrum != null) try makeSpectrum(T, &spectrum, acf, allocator);

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, acf, spectrum, wavefunction); return output;
}

/// Returns the propagators for the k-space grid.
pub fn kgridPropagators(comptime T: type, nstate: u32, kvec: Matrix(T), time_step: T, mass: T, imaginary: bool, allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
    const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var K = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, kvec.rows);

    for (0..kvec.rows) |i| {

        try K.append(try Matrix(Complex(T)).init(nstate, nstate, allocator));

        for (0..K.items[i].rows) |j| {

            for(0..kvec.cols) |k| K.items[i].ptr(j, j).* = K.items[i].at(j, j).add(Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));

            K.items[i].ptr(j, j).* = std.math.complex.exp(K.items[i].at(j, j).mul(Complex(T).init(-0.5 * time_step / mass, 0)).mul(unit));
        }
    }

    return K;
}

/// Function to transform the autocorrelation function to a spectrum.
pub fn makeSpectrum(comptime T: type, spectrum: *Matrix(T), acf: Matrix(T), allocator: std.mem.Allocator) !void {
    var transform = try Vector(Complex(T)).init(spectrum.rows,    allocator); defer transform.deinit();
    var frequency = try Matrix(T         ).init(spectrum.rows, 1, allocator); defer frequency.deinit();

    mpt.kgrid(T, &frequency, 0, asfloat(T, spectrum.rows) * (acf.at(1, 0) - acf.at(0, 0)), @intCast(spectrum.rows));

    for (0..acf.rows) |i| {
        transform.ptr(acf.rows + i).* = Complex(T).init(acf.at(i, 1), -acf.at(i, 2)).mul(Complex(T).init(std.math.exp(-0.001 * acf.at(i, 0) * acf.at(i, 0)), 0));
        transform.ptr(acf.rows - i).* = Complex(T).init(acf.at(i, 1),  acf.at(i, 2)).mul(Complex(T).init(std.math.exp(-0.001 * acf.at(i, 0) * acf.at(i, 0)), 0));
    }

    try ftr.fft(T, transform.sa(), -1);

    for (0..spectrum.rows) |i| {
        spectrum.ptr(i, 0).* = frequency.at(i, 0); spectrum.ptr(i, 1).* = transform.at(i   ).magnitude();
    }

    for (0..spectrum.rows) |i| for (1..spectrum.rows - i) |j| if (spectrum.at(j - 1, 0) > spectrum.at(j, 0)) for (0..spectrum.cols) |k| {
        std.mem.swap(T, spectrum.ptr(j - 1, k), spectrum.ptr(j, k));
    };
}

/// Returns the potential matrices for each point in the space.
pub fn rgridPotentials(comptime T: type, potential: []const u8, rvec: Matrix(T), allocator: std.mem.Allocator) ![3]std.ArrayList(Matrix(Complex(T))) {
    const nstate = try mpt.states(potential);

    var T1 = try Matrix(T).init(nstate, nstate, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nstate, nstate, allocator); defer T2.deinit();

    var U   = try Matrix(T).init(nstate, nstate, allocator); defer   U.deinit();
    var UA  = try Matrix(T).init(nstate, nstate, allocator); defer  UA.deinit();
    var UC  = try Matrix(T).init(nstate, nstate, allocator); defer  UC.deinit();
    var UCP = try Matrix(T).init(nstate, nstate, allocator); defer UCP.deinit();

    var V  = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VA = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VC = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        @memcpy(UCP.data, UC.data);

        try mpt.eval(T, &U, potential, rvec.row(i).vector()); mat.eigh(T, &UA, &UC, U, &T1, &T2);

        if (i > 0) for (0..UC.cols) |j| {
            var overlap: T = 0; for (0..UC.rows) |k| {overlap += UC.at(k, j) * UCP.at(k, j);} if (overlap < 0) for (0..UC.rows) |k| {UC.ptr(k, j).* *= -1;};
        };

        try V.append(try U.complex()); try VA.append(try UA.complex()); try VC.append(try UC.complex());
    }

    return .{V, VA, VC};
}

/// Returns the propagators for the r-space grid.
pub fn rgridPropagators(comptime T: type, VA: std.ArrayList(Matrix(Complex(T))), VC: @TypeOf(VA), rvec: Matrix(T), time_step: T, imaginary: bool, allocator: std.mem.Allocator) !@TypeOf(VA) {
    const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var T1 = try Matrix(Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer T1.deinit();
    var T2 = try Matrix(Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer T2.deinit();

    var R = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        T1.fill(Complex(T).init(0, 0)); for (0..T1.rows) |j| {
            T1.ptr(j, j).* = std.math.complex.exp(VA.items[i].at(j, j).mul(Complex(T).init(-0.5 * time_step, 0)).mul(unit));
        }

        mat.mm(Complex(T), &T2, VC.items[i], T1); mat.mma(Complex(T), &T1, T2, VC.items[i]); try R.append(try T1.clone());
    }

    return R;
}

/// Prints the iteration information.
pub fn printIteration(comptime T: type, i: u32, Ekin: T, Epot: T, r: Vector(T), p: Vector(T), P: Matrix(T), acfi: Complex(T)) !void {
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
            try std.io.getStdOut().writer().print("{s}{d:8.5}", .{if (j == 0) "" else ", ", P.at(j, j)});
        }

        try std.io.getStdOut().writer().print("]\n", .{});
}

/// Writes the results of the quantum dynamics to the output files.
pub fn writeResults(comptime T: type, opt: QuantumDynamicsOptions(T), pop: Matrix(T), ekin: Matrix(T), epot: Matrix(T), etot: Matrix(T), position: Matrix(T), momentum: Matrix(T), acf: Matrix(T), spectrum: Matrix(T), wavefunction: Matrix(T)) !void {
    if (opt.write.autocorrelation_function) |path| try          acf.write(path);
    if (opt.write.kinetic_energy          ) |path| try         ekin.write(path);
    if (opt.write.momentum                ) |path| try     momentum.write(path);
    if (opt.write.population              ) |path| try          pop.write(path);
    if (opt.write.position                ) |path| try     position.write(path);
    if (opt.write.potential_energy        ) |path| try         epot.write(path);
    if (opt.write.spectrum                ) |path| try     spectrum.write(path);
    if (opt.write.total_energy            ) |path| try         etot.write(path);
    if (opt.write.wavefunction            ) |path| try wavefunction.write(path);
}

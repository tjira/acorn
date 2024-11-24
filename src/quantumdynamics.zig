const std = @import("std"); const Complex = std.math.Complex; const gsl_eigen = @cImport(@cInclude("gsl/gsl_eigen.h")); const gsl_fft = @cImport(@cInclude("gsl/gsl_fft_complex.h"));

const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );
const wfn = @import("wavefunction.zig"    );

const Matrix       = @import("matrix.zig"      ).Matrix      ;
const Vector       = @import("vector.zig"      ).Vector      ;
const Wavefunction = @import("wavefunction.zig").Wavefunction;

const asfloat = @import("helper.zig").asfloat;

pub fn QuantumDynamicsOptions(comptime T: type) type {
    return struct {
        const Grid = struct {
            limits: []const T, points: u32
        };
        const InitialConditions = struct {
            position: []const T, momentum: []const T, state: u32, mass: T
        };
        const LogIntervals = struct {
            iteration: u32
        };
        const Write = struct {
            kinetic_energy: ?[]const u8,
            momentum: ?[]const u8,
            population: ?[]const u8,
            position: ?[]const u8,
            potential_energy: ?[]const u8,
            total_energy: ?[]const u8
        };

        adiabatic: bool,
        imaginary: bool,
        iterations: u32,
        time_step: T,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals, write: Write, potential: []const u8
    };
}

pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var pop      = try Matrix(T).init(opt.iterations, mpt.states(opt.potential), allocator);     defer  pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer momentum.deinit(); momentum.fill(0);

    {
        const GSLFT = gsl_fft.gsl_fft_complex_wavetable_alloc(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential))); defer gsl_fft.gsl_fft_complex_wavetable_free(GSLFT);
        const GSLFW = gsl_fft.gsl_fft_complex_workspace_alloc(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential))); defer gsl_fft.gsl_fft_complex_workspace_free(GSLFW);

        var rvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), mpt.dims(opt.potential), allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), mpt.dims(opt.potential), allocator); defer kvec.deinit();

        var r = try Vector(T).init(mpt.dims(opt.potential), allocator); defer r.deinit(); r.fill(0);
        var p = try Vector(T).init(mpt.dims(opt.potential), allocator); defer p.deinit(); p.fill(0);

        const dr = (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1);

        var W  = try Wavefunction(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), opt.grid.points, allocator); defer  W.deinit();
        var WA = try Wavefunction(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), opt.grid.points, allocator); defer WA.deinit();

        var P = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer P.deinit();

        var T1 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, W.ndim), 1, allocator); defer T1.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        const VS = try rgridPotentials(T, opt.potential, rvec, allocator); var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

        const R = try rgridPropagators(T, VA, VC,   rvec, opt.time_step,                              opt.imaginary, allocator); defer R.deinit();
        const K = try kgridPropagators(T, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, opt.imaginary, allocator); defer K.deinit();

        wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.state); W.normalize(dr);

        try std.io.getStdOut().writer().print("\n{s:6} {s:12} {s:12} {s:12}\n", .{"ITER", "EKIN", "EPOT", "ETOT"});

        for (0..opt.iterations) |i| {

            wfn.propagate(T, &W, R, K, &T1, GSLFT, GSLFW); if (opt.imaginary) W.normalize(dr);

            if (opt.adiabatic) wfn.adiabatize(T, &WA, W, VC);

            const Ekin = wfn.ekin(T, W, kvec, opt.initial_conditions.mass, dr, &T1, GSLFT, GSLFW); const Epot: T = wfn.epot(T, W, V, dr);

            wfn.density(T, &P, if (opt.adiabatic) WA else W, dr); wfn.position(T, &r, W, rvec, dr); wfn.momentum(T, &p, W, kvec, dr, &T1, GSLFT, GSLFW);

            if (opt.write.population       != null) for (0..W.nstate) |j| {pop.ptr(i, j).* = P.at(j, j);};
            if (opt.write.position         != null) for (0..W.ndim) |j| {position.ptr(i, j).* = r.at(j);};
            if (opt.write.momentum         != null) for (0..W.ndim) |j| {momentum.ptr(i, j).* = p.at(j);};
            if (opt.write.kinetic_energy   != null) ekin.ptr(i, 0).* = Ekin                              ;
            if (opt.write.potential_energy != null) epot.ptr(i, 0).* = Epot                              ;
            if (opt.write.total_energy     != null) etot.ptr(i, 0).* = Ekin + Epot                       ;

            if (i == 0 or (i + 1) % opt.log_intervals.iteration == 0) {
                try std.io.getStdOut().writer().print("{d:6} {d:12.6} {d:12.6} {d:12.6}\n", .{i + 1, Ekin, Epot, Ekin + Epot});
            }
        }

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();} for (VA.items) |*e| {e.deinit();} for (VC.items) |*e| {e.deinit();}
    }

    for (0..mpt.states(opt.potential)) |i| {
        try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});
    }

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, allocator);
}

fn kgridPropagators(comptime T: type, nstate: u32, kvec: Matrix(T), time_step: T, mass: T, imaginary: bool, allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
    const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var T1 = try Matrix(Complex(T)).init(nstate, nstate, allocator); defer T1.deinit();

    var K = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, kvec.rows);

    for (0..kvec.rows) |i| {

        T1.fill(Complex(T).init(0, 0));

        for (0..T1.rows) |j| for(0..kvec.cols) |k| {T1.ptr(j, j).* = T1.at(j, j).add(Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));};

        for (0..T1.rows) |j| T1.ptr(j, j).* = std.math.complex.exp(T1.at(j, j).mul(Complex(T).init(-0.5 * time_step / mass, 0)).mul(unit));

        try K.append(try T1.clone());
    }

    return K;
}

fn rgridPotentials(comptime T: type, potential: []const u8, rvec: Matrix(T), allocator: std.mem.Allocator) ![3]std.ArrayList(Matrix(Complex(T))) {
    const GSLEW = gsl_eigen.gsl_eigen_symmv_alloc(6 * mpt.states(potential)); defer gsl_eigen.gsl_eigen_symmv_free(GSLEW);

    var U  = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer  U.deinit();
    var UA = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UA.deinit();
    var UC = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UC.deinit();
    var UT = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UT.deinit();

    var V  = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VA = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VC = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        mpt.eval(T, &U, potential, rvec.rowptr(i).vectorptr()); mat.eigh(T, &UA, &UC, U, &UT, GSLEW);

        try V.append(try U.complex()); try VA.append(try UA.complex()); try VC.append(try UC.complex());
    }

    return .{V, VA, VC};
}

fn rgridPropagators(comptime T: type, VA: std.ArrayList(Matrix(Complex(T))), VC: @TypeOf(VA), rvec: Matrix(T), time_step: T, imaginary: bool, allocator: std.mem.Allocator) !@TypeOf(VA) {
    const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var T3 = try Matrix(Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer T3.deinit();
    var T4 = try Matrix(Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer T4.deinit();

    var R = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        T3.fill(Complex(T).init(0, 0)); for (0..T3.rows) |j| {
            T3.ptr(j, j).* = std.math.complex.exp(VA.items[i].at(j, j).mul(Complex(T).init(-0.5 * time_step, 0)).mul(unit));
        }

        T4.fill(Complex(T).init(0, 0)); for (0..T4.rows) |j| for (0..T4.cols) |k| for (0..T4.rows) |l| {
            T4.ptr(j, k).* = T4.at(j, k).add(VC.items[i].at(j, l).mul(T3.at(l, k)));
        };

        T3.fill(Complex(T).init(0, 0)); for (0..T3.rows) |j| for (0..T3.cols) |k| for (0..T3.rows) |l| {
            T3.ptr(j, k).* = T3.at(j, k).add(T4.at(j, l).mul(VC.items[i].at(k, l)));
        };

        try R.append(try T3.clone());
    }

    return R;
}

fn writeResults(comptime T: type, opt: QuantumDynamicsOptions(T), pop: Matrix(T), ekin: Matrix(T), epot: Matrix(T), etot: Matrix(T), position: Matrix(T), momentum: Matrix(T), allocator: std.mem.Allocator) !void {
    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * asfloat(T, opt.iterations));

    if (opt.write.kinetic_energy) |path| {
        var ekin_t = try Matrix(T).init(opt.iterations, ekin.cols + 1, allocator); mat.hjoin(T, &ekin_t, time, ekin); try ekin_t.write(path); ekin_t.deinit();
    }

    if (opt.write.population) |path| {
        var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); mat.hjoin(T, &pop_t, time, pop); try pop_t.write(path); pop_t.deinit();
    }

    if (opt.write.potential_energy) |path| {
        var epot_t = try Matrix(T).init(opt.iterations, epot.cols + 1, allocator); mat.hjoin(T, &epot_t, time, epot); try epot_t.write(path); epot_t.deinit();
    }

    if (opt.write.total_energy) |path| {
        var etot_t = try Matrix(T).init(opt.iterations, etot.cols + 1, allocator); mat.hjoin(T, &etot_t, time, etot); try etot_t.write(path); etot_t.deinit();
    }

    if (opt.write.position) |path| {
        var position_t = try Matrix(T).init(opt.iterations, position.cols + 1, allocator); mat.hjoin(T, &position_t, time, position); try position_t.write(path); position_t.deinit();
    }

    if (opt.write.momentum) |path| {
        var momentum_t = try Matrix(T).init(opt.iterations, momentum.cols + 1, allocator); mat.hjoin(T, &momentum_t, time, momentum); try momentum_t.write(path); momentum_t.deinit();
    }
}

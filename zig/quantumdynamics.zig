const std = @import("std"); const Complex = std.math.Complex; const gsl = @cImport(@cInclude("gsl/gsl_eigen.h"));

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
            population: ?[]const u8
        };

        adiabatic: bool,
        imaginary: bool,
        iterations: u32,
        time_step: T,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals, write: Write, potential: []const u8
    };
}

pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var pop = try Matrix(T).init(opt.iterations, mpt.states(opt.potential), allocator); defer pop.deinit(); pop.fill(0);

    {
        var rvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), mpt.dims(opt.potential), allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), mpt.dims(opt.potential), allocator); defer kvec.deinit();

        const dr = (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1);

        var W  = try Wavefunction(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), opt.grid.points, allocator); defer  W.deinit();
        var WA = try Wavefunction(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), opt.grid.points, allocator); defer WA.deinit();

        var P = try Matrix(T).init(mpt.states(opt.potential), mpt.states(opt.potential), allocator); defer P.deinit();

        var T1 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, W.ndim), 1, allocator); defer T1.deinit();
        var T2 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, W.ndim), 1, allocator); defer T2.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        const VS = try rgridPotentials(T, opt.potential, rvec, allocator); var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

        const R = try rgridPropagators(T, VA, VC,   rvec, opt.time_step,                              opt.imaginary, allocator); defer R.deinit();
        const K = try kgridPropagators(T, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, opt.imaginary, allocator); defer K.deinit();

        wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.state); W.normalize(dr);

        try std.io.getStdOut().writer().print("\n{s:6} {s:12} {s:12} {s:12}\n", .{"ITER", "EKIN", "EPOT", "ETOT"});

        for (0..opt.iterations) |i| {

            wfn.propagate(T, &W, R, K, &T1, &T2); if (opt.imaginary) W.normalize(dr);

            if (opt.adiabatic) wfn.adiabatize(T, &WA, W, VC);

            const Ekin = wfn.ekin(T, W, kvec, opt.initial_conditions.mass, dr, &T1, &T2); const Epot: T = wfn.epot(T, W, V, dr);

            wfn.density(T, &P, if (opt.adiabatic) WA else W, dr); for (0..W.nstate) |j| pop.ptr(i, j).* = P.at(j, j);

            if (i == 0 or (i + 1) % opt.log_intervals.iteration == 0) {
                try std.io.getStdOut().writer().print("{d:6} {d:12.6} {d:12.6} {d:12.6}\n", .{i + 1, Ekin, Epot, Ekin + Epot});
            }
        }

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();} for (VA.items) |*e| {e.deinit();} for (VC.items) |*e| {e.deinit();}
    }

    for (0..mpt.states(opt.potential)) |i| {
        try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});
    }

    try writeResults(T, opt, pop, allocator);
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
    const GSLW = gsl.gsl_eigen_symmv_alloc(mpt.states(potential)); defer gsl.gsl_eigen_symmv_free(GSLW);

    var U  = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer  U.deinit();
    var UA = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UA.deinit();
    var UC = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UC.deinit();

    var V  = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VA = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VC = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        mpt.eval(T, &U, potential, rvec.rowptr(i).vectorptr()); mat.eigh(T, &UA, &UC, U, GSLW);

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

fn writeResults(comptime T: type, opt: QuantumDynamicsOptions(T), pop: Matrix(T), allocator: std.mem.Allocator) !void {
    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * asfloat(T, opt.iterations));

    if (opt.write.population) |path| {
        var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); mat.hjoin(T, &pop_t, time, pop); try pop_t.write(path); pop_t.deinit();
    }
}

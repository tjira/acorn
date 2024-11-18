const std = @import("std"); const Complex = std.math.Complex;

const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );
const wfn = @import("wavefunction.zig"    );

const Matrix       = @import("matrix.zig"      ).Matrix      ;
const Vector       = @import("vector.zig"      ).Vector      ;
const Wavefunction = @import("wavefunction.zig").Wavefunction;

const asfloat = @import("helper.zig").asfloat;

const stdout = std.io.getStdOut().writer();

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

        imaginary: bool,
        iterations: u32,
        time_step: f64,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals, write: Write, potential: mpt.PotentialType(T),
    };
}

pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    var pop = try Matrix(T).init(opt.iterations, try mpt.states(T, opt.potential), allocator); defer pop.deinit(); pop.fill(0);

    {
        var rvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, try mpt.dims(T, opt.potential)), try mpt.dims(T, opt.potential), allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, try mpt.dims(T, opt.potential)), try mpt.dims(T, opt.potential), allocator); defer kvec.deinit();

        const dr = (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1);

        var W = try Wavefunction(T).init(try mpt.dims(T, opt.potential), try mpt.states(T, opt.potential), opt.grid.points, allocator); defer W.deinit();

        var P = try Matrix(T).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer P.deinit();

        var T1 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, W.ndim), 1, allocator); defer T1.deinit();
        var T2 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, W.ndim), 1, allocator); defer T2.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        var R = try rgridPropagators(T, opt, rvec, allocator); defer R.deinit();
        var K = try kgridPropagators(T, opt, kvec, allocator); defer K.deinit();
        var V = try rgridPotentials (T, opt, rvec, allocator); defer V.deinit();

        wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.state); W.normalize(dr);

        for (0..opt.iterations) |i| {

            wfn.propagate(T, &W, R, K, &T1, &T2); if (opt.imaginary) W.normalize(dr);

            const Ekin = wfn.ekin(T, W, kvec, opt.initial_conditions.mass, dr, &T1, &T2); const Epot: T = wfn.epot(T, W, V, dr);

            wfn.density(T, &P, W, dr); for (0..W.nstate) |j| pop.ptr(i, j).* = P.at(j, j);

            if (i == 0 or (i + 1) % opt.log_intervals.iteration == 0) {
                try stdout.print("{d:6} {d:12.6} {d:12.6} {d:12.6}", .{i + 1, Ekin, Epot, Ekin + Epot});
            }

            if (i == 0 or (i + 1) % opt.log_intervals.iteration == 0) {
                for (0..W.nstate) |j| try stdout.print(" {d:12.6}", .{P.at(j, j)}); try stdout.print("\n", .{});
            }
        }

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();}
    }

    for (0..try mpt.states(T, opt.potential)) |i| {
        try stdout.print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});
    }

    try writeResults(T, opt, pop, allocator);
}

pub fn kgridPropagators(comptime T: type, opt: QuantumDynamicsOptions(T), kvec: Matrix(T), allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
    const unit = Complex(T).init(if (opt.imaginary) 1 else 0, if (opt.imaginary) 0 else 1);

    var T1 = try Matrix(Complex(T)).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer T1.deinit();

    var K = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, kvec.rows);

    for (0..kvec.rows) |i| {

        T1.fill(Complex(T).init(0, 0));

        for (0..T1.rows) |j| for(0..kvec.cols) |k| {T1.ptr(j, j).* = T1.at(j, j).add(Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));};

        for (0..T1.rows) |j| T1.ptr(j, j).* = std.math.complex.exp(T1.at(j, j).mul(Complex(T).init(-0.5 * opt.time_step / opt.initial_conditions.mass, 0)).mul(unit));

        try K.append(try T1.clone());
    }

    return K;
}

pub fn rgridPotentials(comptime T: type, opt: QuantumDynamicsOptions(T), rvec: Matrix(T), allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
    var U  = try Matrix(T         ).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer  U.deinit();
    var T1 = try Matrix(Complex(T)).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer T1.deinit();
    var r  = try Vector(T         ).init(try mpt.dims  (T, opt.potential),                                   allocator); defer  r.deinit();

    var V = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {
        
        for (0..r.rows) |j| {r.ptr(j).* = rvec.at(i, j);} opt.potential(T, &U, r);

        for (0..U.rows) |j| for (0..U.rows) |k| {
            T1.ptr(j, k).* = Complex(T).init(U.at(j, k), 0);
        };

        try V.append(try T1.clone());
    }

    return V;
}

pub fn rgridPropagators(comptime T: type, opt: QuantumDynamicsOptions(T), rvec: Matrix(T), allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
    const unit = Complex(T).init(if (opt.imaginary) 1 else 0, if (opt.imaginary) 0 else 1);

    var U  = try Matrix(T         ).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer  U.deinit();
    var UA = try Matrix(T         ).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer UA.deinit();
    var UC = try Matrix(T         ).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer UC.deinit();
    var T1 = try Matrix(T         ).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer T1.deinit();
    var T2 = try Matrix(T         ).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer T2.deinit();
    var T3 = try Matrix(Complex(T)).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer T3.deinit();
    var T4 = try Matrix(Complex(T)).init(try mpt.states(T, opt.potential), try mpt.states(T, opt.potential), allocator); defer T4.deinit();
    var r  = try Vector(T         ).init(try mpt.dims  (T, opt.potential),                                   allocator); defer  r.deinit();

    var R = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, opt.grid.points);

    for (0..rvec.rows) |i| {

        for (0..r.rows) |j| {r.ptr(j).* = rvec.at(i, j);} opt.potential(T, &U, r); mat.eigh(T, &UA, &UC, U, 1e-12, &T1, &T2); T3.fill(Complex(T).init(0, 0));

        for (0..U.rows) |j| T3.ptr(j, j).* = std.math.complex.exp(Complex(T).init(UA.at(j, j), 0).mul(Complex(T).init(-0.5 * opt.time_step, 0)).mul(unit));

        T4.fill(Complex(T).init(0, 0)); for (0..U.rows) |j| for (0..U.cols) |k| for (0..U.rows) |l| {
            T4.ptr(j, k).* = T4.at(j, k).add(Complex(T).init(UC.at(j, l), 0).mul(T3.at(l, k)));
        };

        T3.fill(Complex(T).init(0, 0)); for (0..U.rows) |j| for (0..U.cols) |k| for (0..U.rows) |l| {
            T3.ptr(j, k).* = T3.at(j, k).add(T4.at(j, l).mul(Complex(T).init(UC.at(k, l), 0)));
        };

        try R.append(try T3.clone());
    }

    return R;
}

fn writeResults(comptime T: type, opt: QuantumDynamicsOptions(T), pop: Matrix(T), allocator: std.mem.Allocator) !void {
    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * opt.iterations);

    if (opt.write.population) |path| {
        var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); time.hjoin(&pop_t, pop); try pop_t.write(path); pop_t.deinit();
    }
}

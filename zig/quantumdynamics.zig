const std = @import("std");

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

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

        imaginary: bool,
        iterations: u32,
        time_step: f64,

        grid: Grid, initial_conditions: InitialConditions, log_intervals: LogIntervals, potential: mpt.PotentialType(T),
    };
}

pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), allocator: std.mem.Allocator) !void {
    const ndim = try mpt.dims(T, opt.potential); const nstate = try mpt.states(T, opt.potential);

    var pop = try Matrix(T).init(opt.iterations, nstate, allocator); defer pop.deinit(); pop.fill(0);

    {
        var rvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, ndim), ndim, allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, ndim), ndim, allocator); defer kvec.deinit();

        var P = try Matrix(T).init(nstate, nstate, allocator); defer P.deinit();

        var W  = try Matrix(std.math.Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), nstate, allocator); defer  W.deinit();
        var T1 = try Matrix(std.math.Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), nstate, allocator); defer T1.deinit();
        var T2 = try Matrix(std.math.Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), nstate, allocator); defer T2.deinit();

        var T3 = try Matrix(std.math.Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), 1, allocator); defer T3.deinit();
        var T4 = try Matrix(std.math.Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), 1, allocator); defer T4.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        var R = try rgridPropagators(T, opt, rvec, allocator); defer R.deinit();
        var K = try kgridPropagators(T, opt, kvec, allocator); defer K.deinit();
        var V = try rgridPotentials (T, opt, rvec, allocator); defer V.deinit();

        wfnGuess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.state); const dr = rvec.at(1, ndim - 1) - rvec.at(0, ndim - 1);

        // for (0..R.items.len) |i| try stdout.print("{d:12.6} {d:12.6}\n", .{K.items[i].at(0, 0).re, K.items[i].at(0, 0).im});

        for (0..opt.iterations) |i| {

            wfnPropagate(T, &W, R, K, &T3, &T4);

            // print the wfn
            // for (0..W.rows) |j| try stdout.print("{d:12.6} {d:12.6} {d:12.6} {d:12.6}\n", .{W.at(j, 0).re, W.at(j, 0).im, W.at(j, 1).re, W.at(j, 1).im});

            if (opt.imaginary) {const norm = wfnNorm(T, W, dr); for (0..W.rows) |j| for (0..W.cols) |k| {W.ptr(j, k).* = W.at(j, k).div(std.math.Complex(T).init(norm, 0));};}

            const Ekin = wfnKineticEnergy(T, W, kvec, opt.initial_conditions.mass, dr, &T3, &T4); const Epot: T = wfnPotentialEnergy(T, W, V, dr);

            wfnDensity(T, &P, W, dr); for (0..nstate) |j| pop.ptr(i, j).* = P.at(j, j);

            try stdout.print("{d:6} {d:12.6} {d:12.6} {d:12.6}", .{i + 1, Ekin, Epot, Ekin + Epot});

            for (0..nstate) |j| try stdout.print(" {d:12.6}", .{P.at(j, j)}); try stdout.print("\n", .{});
        }

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();}
    }

    try writeResults(T, opt, pop, allocator);
}

pub fn kgridPropagators(comptime T: type, opt: QuantumDynamicsOptions(T), kvec: Matrix(T), allocator: std.mem.Allocator) !std.ArrayList(Matrix(std.math.Complex(T))) {
    const ndim = try mpt.dims(T, opt.potential); const nstate = try mpt.states(T, opt.potential);

    const unit = std.math.Complex(T).init(if (opt.imaginary) 1 else 0, if (opt.imaginary) 0 else 1);

    var T1 = try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator); defer T1.deinit();

    var K = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, opt.grid.points);

    for (0..opt.grid.points) |i| {

        T1.fill(std.math.Complex(T).init(0, 0));

        for (0..nstate) |j| for(0..ndim) |k| {T1.ptr(j, j).* = T1.at(j, j).add(std.math.Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));};

        for (0..nstate) |j| T1.ptr(j, j).* = std.math.complex.exp(T1.at(j, j).mul(std.math.Complex(T).init(-0.5 * opt.time_step / opt.initial_conditions.mass, 0)).mul(unit));

        try K.append(try T1.clone());
    }

    return K;
}

pub fn rgridPotentials(comptime T: type, opt: QuantumDynamicsOptions(T), rvec: Matrix(T), allocator: std.mem.Allocator) !std.ArrayList(Matrix(std.math.Complex(T))) {
    const ndim = try mpt.dims(T, opt.potential); const nstate = try mpt.states(T, opt.potential); var r = try Vector(T).init(ndim, allocator); defer r.deinit();

    var U  = try Matrix(T).init(nstate, nstate, allocator); defer  U.deinit();

    var T1 = try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator); defer T1.deinit();

    var V = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, opt.grid.points);

    for (0..opt.grid.points) |i| {
        
        for (0..ndim) |j| {r.ptr(j).* = rvec.at(i, j);} opt.potential(T, &U, r);

        for (0..nstate) |j| for (0..nstate) |k| {T1.ptr(j, k).* = std.math.Complex(T).init(U.at(j, k), 0);};

        try V.append(try T1.clone());
    }

    return V;
}

pub fn rgridPropagators(comptime T: type, opt: QuantumDynamicsOptions(T), rvec: Matrix(T), allocator: std.mem.Allocator) !std.ArrayList(Matrix(std.math.Complex(T))) {
    const ndim = try mpt.dims(T, opt.potential); const nstate = try mpt.states(T, opt.potential); var r = try Vector(T).init(ndim, allocator); defer r.deinit();

    const unit = std.math.Complex(T).init(if (opt.imaginary) 1 else 0, if (opt.imaginary) 0 else 1);

    var U  = try Matrix(T).init(nstate, nstate, allocator); defer  U.deinit();
    var UA = try Matrix(T).init(nstate, nstate, allocator); defer UA.deinit();
    var UC = try Matrix(T).init(nstate, nstate, allocator); defer UC.deinit();
    var T1 = try Matrix(T).init(nstate, nstate, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nstate, nstate, allocator); defer T2.deinit();

    var T3 = try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator); defer T3.deinit();
    var T4 = try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator); defer T4.deinit();
    var T5 = try Matrix(std.math.Complex(T)).init(nstate, nstate, allocator); defer T5.deinit();

    var R = try std.ArrayList(Matrix(std.math.Complex(T))).initCapacity(allocator, opt.grid.points);

    for (0..opt.grid.points) |i| {

        for (0..ndim) |j| {r.ptr(j).* = rvec.at(i, j);} opt.potential(T, &U, r); mat.eigh(T, &UA, &UC, U, 1e-12, &T1, &T2);

        for (0..nstate) |j| for (0..nstate) |k| {T3.ptr(j, k).* = std.math.Complex(T).init(UA.at(j, k), 0);};
        for (0..nstate) |j| for (0..nstate) |k| {T4.ptr(j, k).* = std.math.Complex(T).init(UC.at(j, k), 0);};

        for (0..nstate) |j| T3.ptr(j, j).* = std.math.complex.exp(T3.at(j, j).mul(std.math.Complex(T).init(-0.5 * opt.time_step, 0)).mul(unit));

        mat.cmm(std.math.Complex(T), &T5, T4, T3); mat.transpose(std.math.Complex(T), &T3, T4); mat.cmm(std.math.Complex(T), &T4, T5, T3);

        try R.append(try T4.clone());
    }

    return R;
}

pub fn wfnDensity(comptime T: type, P: *Matrix(T), W: Matrix(std.math.Complex(T)), dr: T) void {
    P.fill(0); var pij = std.math.Complex(T).init(0, 0);

    for (0..W.cols) |i| for (0..W.cols) |j| {
        for (0..W.rows) |k| {pij = pij.add(W.at(k, i).mul(std.math.complex.conj(W.at(k, j))));} P.ptr(i, j).* = std.math.complex.abs(pij) * dr;
    };
}

pub fn wfnGuess(comptime T: type, W: *Matrix(std.math.Complex(T)), rvec: Matrix(T), position: []const T, momentum: []const T, state: u32) void {
    W.fill(std.math.Complex(T).init(0, 0));

    for (0..W.rows) |i| for (0..W.cols) |j| if (j == state) {
        W.ptr(i, j).* = std.math.Complex(T).init(std.math.exp(-(rvec.at(i, 0) - position[0]) * (rvec.at(i, 0) - position[0])), 0);
    };

    for (0..W.rows) |i| for (0..W.cols) |j| if (j == state) {
        W.ptr(i, j).* = W.at(i, j).mul(std.math.complex.exp(std.math.Complex(T).init(0, rvec.at(i, 0) * momentum[0])));
    };

    const norm = wfnNorm(T, W.*, rvec.at(1, 0) - rvec.at(0, 0)); for (0..W.rows) |i| for (0..W.cols) |j| {W.ptr(i, j).* = W.at(i, j).div(std.math.Complex(T).init(norm, 0));};
}

pub fn wfnKineticEnergy(comptime T: type, W: Matrix(std.math.Complex(T)), kvec: Matrix(T), mass: T, dr: T, T1: *Matrix(std.math.Complex(T)), T2: @TypeOf(T1)) T {
    var Ekin: T = 0;

    for (0..W.cols) |i| {

        W.col(T1, i); ftr.dft(T, T2.data, T1.data, &[_]usize{W.rows}, -1);

        for (0..W.rows) |j| T2.ptr(j, 0).* = T2.at(j, 0).mul(std.math.Complex(T).init(kvec.at(j, 0) * kvec.at(j, 0), 0));

        ftr.dft(T, T1.data, T2.data, &[_]usize{W.rows},  1);

        for (0..W.rows) |j| Ekin += (T1.at(j, 0).mul(std.math.complex.conj(W.at(j, i)))).re;
    }

    return 0.5 * Ekin / mass * dr;
}

pub fn wfnNorm(comptime T: type, W: Matrix(std.math.Complex(T)), dr: T) T {
    var norm = std.math.Complex(T).init(0, 0);

    for (0..W.cols) |i| for (0..W.cols) |j| for (0..W.rows) |k| {
        norm = norm.add(std.math.complex.conj(W.at(k, i)).mul(W.at(k, j)));
    };

    return std.math.sqrt(std.math.complex.abs(norm) * dr);
}

pub fn wfnPotentialEnergy(comptime T: type, W: Matrix(std.math.Complex(T)), V: std.ArrayList(Matrix(std.math.Complex(T))), dr: T) T {
    var Epot: T = 0;

    for (0..W.cols) |i| for (0..W.cols) |j| for (0..W.rows) |k| {
        Epot += std.math.complex.abs(std.math.complex.conj(W.at(k, i)).mul(V.items[k].at(i, j).mul(W.at(k, j))));
    };

    return Epot * dr;
}

pub fn wfnPropagate(comptime T: type, W: *Matrix(std.math.Complex(T)), R: std.ArrayList(Matrix(std.math.Complex(T))), K: @TypeOf(R), T1: @TypeOf(W), T2: @TypeOf(W)) void {
    var value = std.math.Complex(T).init(0, 0);

    // for (0..W.rows) |j| W.ptr(j, 0).* = W.at(j, 0).mul(R.items[j].at(0, 0));
    for (0..W.rows) |i| {
        for (0..W.cols) |j| {
            value = std.math.Complex(T).init(0, 0);
            for (0..W.cols) |k| {
                value = value.add(W.at(i, k).mul(R.items[i].at(j, k)));
            }
            W.ptr(i, j).* = value;
        }
    }

    // _ = T1; _ = T2; _ = K;

    // ftr.dft(T, T1.data, W.data, &[_]usize{W.rows}, -1);
    for (0..W.cols) |j| {W.col(T1, j); ftr.dft(T, T2.data, T1.data, &[_]usize{W.rows}, -1); W.setcol(j, T2.*);}

    // for (0..W.rows) |j| T1.ptr(j, 0).* = T1.at(j, 0).mul(K.items[j].at(0, 0));
    for (0..W.rows) |i| {
        for (0..W.cols) |j| {
            value = std.math.Complex(T).init(0, 0);
            for (0..W.cols) |k| {
                value = value.add(W.at(i, k).mul(K.items[i].at(j, k)));
            }
            W.ptr(i, j).* = value;
        }
    }

    // ftr.dft(T, W.data, T1.data, &[_]usize{W.rows},  1);
    for (0..W.cols) |j| {W.col(T1, j); ftr.dft(T, T2.data, T1.data, &[_]usize{W.rows},  1); W.setcol(j, T2.*);}

    // for (0..W.rows) |j| W.ptr(j, 0).* = W.at(j, 0).mul(R.items[j].at(0, 0));
    for (0..W.rows) |i| {
        for (0..W.cols) |j| {
            value = std.math.Complex(T).init(0, 0);
            for (0..W.cols) |k| {
                value = value.add(W.at(i, k).mul(R.items[i].at(j, k)));
            }
            W.ptr(i, j).* = value;
        }
    }
}

fn writeResults(comptime T: type, opt: QuantumDynamicsOptions(T), pop: Matrix(T), allocator: std.mem.Allocator) !void {
    const time = try Matrix(T).init(opt.iterations, 1, allocator); defer time.deinit(); time.linspace(opt.time_step, opt.time_step * opt.iterations);

    var pop_t = try Matrix(T).init(opt.iterations, pop.cols + 1, allocator); time.hjoin(&pop_t, pop); try pop_t.write("POPULATION_EXACT.mat"); pop_t.deinit();
}

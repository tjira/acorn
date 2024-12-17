const std = @import("std"); const Complex = std.math.Complex;

const mat = @import("matrix.zig"        );
const mpt = @import("modelpotential.zig");
const vec = @import("vector.zig"        );
const wfn = @import("wavefunction.zig"  );

const Matrix       = @import("matrix.zig"      ).Matrix      ;
const Vector       = @import("vector.zig"      ).Vector      ;
const Wavefunction = @import("wavefunction.zig").Wavefunction;

const asfloat = @import("helper.zig").asfloat;

pub fn QuantumDynamicsOptions(comptime T: type) type {
    return struct {
        const Grid = struct {
            limits: []const T = &[_]f64{-16, 32}, points: u32 = 512
        };
        const InitialConditions = struct {
            position: []const T = &[_]f64{-10}, momentum: []const T = &[_]f64{15}, gamma: T = 2, state: u32 = 1, mass: T = 2000
        };
        const LogIntervals = struct {
            iteration: u32 = 1
        };
        const Write = struct {
            kinetic_energy: ?[]const u8 = null,
            momentum: ?[]const u8 = null,
            population: ?[]const u8 = null,
            position: ?[]const u8 = null,
            potential_energy: ?[]const u8 = null,
            total_energy: ?[]const u8 = null
        };

        adiabatic: bool = true,
        imaginary: bool = false,
        iterations: u32 = 300,
        time_step: T = 10,

        grid: Grid = .{}, initial_conditions: InitialConditions = .{}, log_intervals: LogIntervals = .{}, write: Write = .{}, potential: []const u8 = "tully1D_1"
    };
}

pub fn QuantumDynamicsOutput(comptime T: type) type {
    return struct {
        P: Matrix(T),
        r: Vector(T),
        p: Vector(T),
        Ekin: T,
        Epot: T,

        pub fn init(ndim: usize, nstate: usize, allocator: std.mem.Allocator) !QuantumDynamicsOutput(T) {
            return QuantumDynamicsOutput(T){
                .P = try Matrix(T).init(nstate, nstate, allocator),
                .r = try Vector(T).init(ndim, allocator),
                .p = try Vector(T).init(ndim, allocator),
                .Ekin = undefined,
                .Epot = undefined
            };
        }
        pub fn deinit(self: QuantumDynamicsOutput(T)) void {
            self.P.deinit(); self.r.deinit(); self.p.deinit();
        }
    };
}

pub fn run(comptime T: type, opt: QuantumDynamicsOptions(T), print: bool, allocator: std.mem.Allocator) !QuantumDynamicsOutput(T) {
    const ndim = try mpt.dims(opt.potential); const nstate = try mpt.states(opt.potential);

    var output = try QuantumDynamicsOutput(T).init(ndim, nstate, allocator);

    var pop      = try Matrix(T).init(opt.iterations, nstate, allocator);     defer  pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations, ndim  , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations, ndim  , allocator); defer momentum.deinit(); momentum.fill(0);

    {
        var T1 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), 1, allocator); defer T1.deinit();
        var T2 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, ndim), 1, allocator); defer T2.deinit();

        var rvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, ndim), ndim, allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, ndim), ndim, allocator); defer kvec.deinit();

        var r = try Vector(T).init(ndim, allocator); defer r.deinit(); r.fill(0);
        var p = try Vector(T).init(ndim, allocator); defer p.deinit(); p.fill(0);

        const dr = (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1);

        var W  = try Wavefunction(T).init(ndim, nstate, opt.grid.points, allocator); defer  W.deinit();
        var WA = try Wavefunction(T).init(ndim, nstate, opt.grid.points, allocator); defer WA.deinit();

        var P  = try Matrix(T).init(nstate, nstate, allocator); defer  P.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        const VS = try rgridPotentials(T, opt.potential, rvec, allocator); var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

        const R = try rgridPropagators(T, VA, VC,   rvec, opt.time_step,                              opt.imaginary, allocator); defer R.deinit();
        const K = try kgridPropagators(T, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, opt.imaginary, allocator); defer K.deinit();

        wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.gamma, opt.initial_conditions.state); wfn.normalize(T, &W, dr);

        if (print) try std.io.getStdOut().writer().print("\n{s:6} {s:12} {s:12} {s:12}", .{"ITER", "ET1N", "EPOT", "ETOT"});

        if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"POSITION"  });}
        if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"MOMENTUM"  });}
        if (print) {if (W.nstate > 1) for (0..W.nstate - 1) |_| {try std.io.getStdOut().writer().print(" " ** 10, .{});}; try std.io.getStdOut().writer().print(" {s:10}\n", .{"POPULATION"});}

        for (0..opt.iterations) |i| {

            wfn.propagate(T, &W, R, K, &T1, &T2); if (opt.imaginary) wfn.normalize(T, &W, dr);

            if (opt.adiabatic) wfn.adiabatize(T, &WA, W, VC);

            const Ekin = wfn.ekin(T, W, kvec, opt.initial_conditions.mass, dr, &T1, &T2); const Epot: T = wfn.epot(T, W, V, dr);

            wfn.density(T, &P, if (opt.adiabatic) WA else W, dr); wfn.position(T, &r, W, rvec, dr); wfn.momentum(T, &p, W, kvec, dr, &T1, &T2);

            if (opt.write.population       != null) for (0..W.nstate) |j| {pop.ptr(i, j).* = P.at(j, j);};
            if (opt.write.position         != null) for (0..W.ndim) |j| {position.ptr(i, j).* = r.at(j);};
            if (opt.write.momentum         != null) for (0..W.ndim) |j| {momentum.ptr(i, j).* = p.at(j);};
            if (opt.write.kinetic_energy   != null) ekin.ptr(i, 0).* = Ekin                              ;
            if (opt.write.potential_energy != null) epot.ptr(i, 0).* = Epot                              ;
            if (opt.write.total_energy     != null) etot.ptr(i, 0).* = Ekin + Epot                       ;

            if (i == opt.iterations - 1) {
                @memcpy(output.P.data, P.data); @memcpy(output.r.data, r.data); @memcpy(output.p.data, p.data); output.Ekin = Ekin; output.Epot = Epot;
            }

            if (print and (i == 0 or (i + 1) % opt.log_intervals.iteration == 0)) try printIteration(T, @intCast(i), Ekin, Epot, r, p, P);
        }

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();} for (VA.items) |*e| {e.deinit();} for (VC.items) |*e| {e.deinit();}
    }

    for (0..nstate) |i| {
        if (print) {try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});}
    }

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, allocator); return output;
}

fn kgridPropagators(comptime T: type, nstate: u32, kvec: Matrix(T), time_step: T, mass: T, imaginary: bool, allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
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

fn rgridPotentials(comptime T: type, potential: []const u8, rvec: Matrix(T), allocator: std.mem.Allocator) ![3]std.ArrayList(Matrix(Complex(T))) {
    const nstate = try mpt.states(potential);

    var T1 = try Matrix(T).init(nstate, nstate, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nstate, nstate, allocator); defer T2.deinit();

    var U  = try Matrix(T).init(nstate, nstate, allocator); defer  U.deinit();
    var UA = try Matrix(T).init(nstate, nstate, allocator); defer UA.deinit();
    var UC = try Matrix(T).init(nstate, nstate, allocator); defer UC.deinit();

    var V  = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VA = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VC = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        try mpt.eval(T, &U, potential, rvec.row(i).vector()); mat.eigh(T, &UA, &UC, U, &T1, &T2);

        try V.append(try U.complex()); try VA.append(try UA.complex()); try VC.append(try UC.complex());
    }

    return .{V, VA, VC};
}

fn rgridPropagators(comptime T: type, VA: std.ArrayList(Matrix(Complex(T))), VC: @TypeOf(VA), rvec: Matrix(T), time_step: T, imaginary: bool, allocator: std.mem.Allocator) !@TypeOf(VA) {
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

fn printIteration(comptime T: type, i: u32, Ekin: T, Epot: T, r: Vector(T), p: Vector(T), P: Matrix(T)) !void {
        try std.io.getStdOut().writer().print("{d:6} {d:12.6} {d:12.6} {d:12.6} [", .{i + 1, Ekin, Epot, Ekin + Epot});

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

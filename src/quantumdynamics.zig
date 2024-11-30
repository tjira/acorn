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
    var output = try QuantumDynamicsOutput(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), allocator);

    var pop      = try Matrix(T).init(opt.iterations, mpt.states(opt.potential), allocator);     defer  pop.deinit();      pop.fill(0);
    var ekin     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer ekin.deinit();     ekin.fill(0);
    var epot     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer epot.deinit();     epot.fill(0);
    var etot     = try Matrix(T).init(opt.iterations, 1                        , allocator);     defer etot.deinit();     etot.fill(0);
    var position = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer position.deinit(); position.fill(0);
    var momentum = try Matrix(T).init(opt.iterations, mpt.dims(opt.potential)  , allocator); defer momentum.deinit(); momentum.fill(0);

    {
        var T1 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), 1, allocator); defer T1.deinit();
        var T2 = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), 1, allocator); defer T2.deinit();

        var rvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), mpt.dims(opt.potential), allocator); defer rvec.deinit();
        var kvec = try Matrix(T).init(std.math.pow(u32, opt.grid.points, mpt.dims(opt.potential)), mpt.dims(opt.potential), allocator); defer kvec.deinit();

        var r = try Vector(T).init(mpt.dims(opt.potential), allocator); defer r.deinit(); r.fill(0);
        var p = try Vector(T).init(mpt.dims(opt.potential), allocator); defer p.deinit(); p.fill(0);

        const dr = (opt.grid.limits[1] - opt.grid.limits[0]) / asfloat(T, opt.grid.points - 1);

        var W  = try Wavefunction(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), opt.grid.points, allocator); defer  W.deinit();
        var WA = try Wavefunction(T).init(mpt.dims(opt.potential), mpt.states(opt.potential), opt.grid.points, allocator); defer WA.deinit();

        var WS = try Matrix(Complex(T)).init(std.math.pow(u32, opt.grid.points, W.ndim), 1,                         allocator); defer WS.deinit();
        var P  = try Matrix(T         ).init(mpt.states(opt.potential),                  mpt.states(opt.potential), allocator); defer  P.deinit();

        mpt.rgrid(T, &rvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);
        mpt.kgrid(T, &kvec, opt.grid.limits[0], opt.grid.limits[1], opt.grid.points);

        const VS = try rgridPotentials(T, opt.potential, rvec, allocator); var V = VS[0]; var VA = VS[1]; var VC = VS[2]; defer V.deinit(); defer VA.deinit(); defer VC.deinit();

        const R = try rgridPropagators(T, VA, VC,   rvec, opt.time_step,                              opt.imaginary, allocator); defer R.deinit();
        const K = try kgridPropagators(T, W.nstate, kvec, opt.time_step, opt.initial_conditions.mass, opt.imaginary, allocator); defer K.deinit();

        wfn.guess(T, &W, rvec, opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.state); W.normalize(dr);

        if (print) try std.io.getStdOut().writer().print("\n{s:6} {s:12} {s:12} {s:12}", .{"ITER", "EKIN", "EPOT", "ETOT"});

        if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"POSITION"  });}
        if (print) {if (W.ndim   > 1) for (0..W.ndim   - 1) |_| {try std.io.getStdOut().writer().print(" " ** 11, .{});}; try std.io.getStdOut().writer().print(" {s:11}",   .{"MOMENTUM"  });}
        if (print) {if (W.nstate > 1) for (0..W.nstate - 1) |_| {try std.io.getStdOut().writer().print(" " ** 10, .{});}; try std.io.getStdOut().writer().print(" {s:10}\n", .{"POPULATION"});}

        for (0..opt.iterations) |i| {

            wfn.propagate(T, &W, R, K, &T1, &T2); if (opt.imaginary) W.normalize(dr);

            if (opt.adiabatic) wfn.adiabatize(T, &WA, W, VC);

            const Ekin = wfn.ekin(T, W, kvec, opt.initial_conditions.mass, dr, &T1, &T2); const Epot: T = wfn.epot(T, W, V, dr);

            wfn.density(T, &P, if (opt.adiabatic) WA else W, dr); wfn.position(T, &r, W, rvec, dr); wfn.momentum(T, &p, W, kvec, dr, &T1, &T2);

            if (opt.write.population       != null) for (0..W.nstate) |j| {pop.ptr(i, j).* = P.at(j, j);};
            if (opt.write.position         != null) for (0..W.ndim) |j| {position.ptr(i, j).* = r.at(j);};
            if (opt.write.momentum         != null) for (0..W.ndim) |j| {momentum.ptr(i, j).* = p.at(j);};
            if (opt.write.kinetic_energy   != null) ekin.ptr(i, 0).* = Ekin                              ;
            if (opt.write.potential_energy != null) epot.ptr(i, 0).* = Epot                              ;
            if (opt.write.total_energy     != null) etot.ptr(i, 0).* = Ekin + Epot                       ;

            @memcpy(output.P.data, P.data); @memcpy(output.r.data, r.data); @memcpy(output.p.data, p.data); output.Ekin = Ekin; output.Epot = Epot;

            if (print and (i == 0 or (i + 1) % opt.log_intervals.iteration == 0)) try printIteration(T, @intCast(i), Ekin, Epot, r, p, P);
        }

        for (R.items) |*e| {e.deinit();} for (K.items) |*e| {e.deinit();} for (V.items) |*e| {e.deinit();} for (VA.items) |*e| {e.deinit();} for (VC.items) |*e| {e.deinit();}
    }

    for (0..mpt.states(opt.potential)) |i| {
        if (print) {try std.io.getStdOut().writer().print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6}\n", .{if (i == 0) "\n" else "", i, pop.at(opt.iterations - 1, i)});}
    }

    try writeResults(T, opt, pop, ekin, epot, etot, position, momentum, allocator); return output;
}

fn kgridPropagators(comptime T: type, nstate: u32, kvec: Matrix(T), time_step: T, mass: T, imaginary: bool, allocator: std.mem.Allocator) !std.ArrayList(Matrix(Complex(T))) {
    const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var KI = try Matrix(Complex(T)).init(nstate, nstate, allocator); defer KI.deinit();

    var K = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, kvec.rows);

    for (0..kvec.rows) |i| {

        KI.fill(Complex(T).init(0, 0)); for (0..KI.rows) |j| {

            for(0..kvec.cols) |k| KI.ptr(j, j).* = KI.at(j, j).add(Complex(T).init(kvec.at(i, k) * kvec.at(i, k), 0));

            KI.ptr(j, j).* = std.math.complex.exp(KI.at(j, j).mul(Complex(T).init(-0.5 * time_step / mass, 0)).mul(unit));
        }

        try K.append(try KI.clone());
    }

    return K;
}

fn rgridPotentials(comptime T: type, potential: []const u8, rvec: Matrix(T), allocator: std.mem.Allocator) ![3]std.ArrayList(Matrix(Complex(T))) {
    var T1 = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer T2.deinit();

    var U  = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer  U.deinit();
    var UA = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UA.deinit();
    var UC = try Matrix(T).init(mpt.states(potential), mpt.states(potential), allocator); defer UC.deinit();

    var V  = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VA = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);
    var VC = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        mpt.eval(T, &U, potential, rvec.rowptr(i).vectorptr()); mat.eigh(T, &UA, &UC, U, &T1, &T2);

        try V.append(try U.complex()); try VA.append(try UA.complex()); try VC.append(try UC.complex());
    }

    return .{V, VA, VC};
}

fn rgridPropagators(comptime T: type, VA: std.ArrayList(Matrix(Complex(T))), VC: @TypeOf(VA), rvec: Matrix(T), time_step: T, imaginary: bool, allocator: std.mem.Allocator) !@TypeOf(VA) {
    const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

    var RI1 = try Matrix(Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer RI1.deinit();
    var RI2 = try Matrix(Complex(T)).init(VA.items[0].rows, VA.items[0].rows, allocator); defer RI2.deinit();

    var R = try std.ArrayList(Matrix(Complex(T))).initCapacity(allocator, rvec.rows);

    for (0..rvec.rows) |i| {

        RI1.fill(Complex(T).init(0, 0)); for (0..RI1.rows) |j| {
            RI1.ptr(j, j).* = std.math.complex.exp(VA.items[i].at(j, j).mul(Complex(T).init(-0.5 * time_step, 0)).mul(unit));
        }

        mat.cmm(T, &RI2, VC.items[i], RI1); mat.cmma(T, &RI1, RI2, VC.items[i]); try R.append(try RI1.clone());
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

// TESTS ===============================================================================================================================================================================================

test "doubleState1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.68939163804408; P.ptr(0, 1).* = 0.00018765157757;
    P.ptr(1, 0).* = 0.00018765157757; P.ptr(1, 1).* = 0.31060836195428;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "doubleState1D_1"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "doubleState1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.96070275892968; P.ptr(0, 1).* = 0.01414658057574;
    P.ptr(1, 0).* = 0.01414658057574; P.ptr(1, 1).* = 0.03929724106854;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "doubleState1D_2"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "tully1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.41050635034636; P.ptr(0, 1).* = 0.03954870766906;
    P.ptr(1, 0).* = 0.03954870766906; P.ptr(1, 1).* = 0.58949364965213;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tully1D_1"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "tully1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.02483440861074; P.ptr(0, 1).* = 0.00092889832192;
    P.ptr(1, 0).* = 0.00092889832192; P.ptr(1, 1).* = 0.97516559138780;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tully1D_2"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "tully1D_3" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.76438466406603; P.ptr(0, 1).* = 0.17171831933474;
    P.ptr(1, 0).* = 0.17171831933474; P.ptr(1, 1).* = 0.23561533593184;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tully1D_3"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "tripleState1D_1" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.09406137151353; P.ptr(0, 1).* = 0.02243393031679; P.ptr(0, 2).* = 0.00177590480664;
    P.ptr(1, 0).* = 0.02243393031679; P.ptr(1, 1).* = 0.09966537676404; P.ptr(1, 2).* = 0.03403002401355;
    P.ptr(2, 0).* = 0.00177590480664; P.ptr(2, 1).* = 0.03403002401355; P.ptr(2, 2).* = 0.80627325172103;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tripleState1D_1"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "tripleState1D_2" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.80243373396786; P.ptr(0, 1).* = 0.09151956326301; P.ptr(0, 2).* = 0.02478818099951;
    P.ptr(1, 0).* = 0.09151956326301; P.ptr(1, 1).* = 0.04418760360549; P.ptr(1, 2).* = 0.04060912220532;
    P.ptr(2, 0).* = 0.02478818099951; P.ptr(2, 1).* = 0.04060912220532; P.ptr(2, 2).* = 0.15337866242496;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tripleState1D_2"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "tripleState1D_3" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.72007774176220; P.ptr(0, 1).* = 0.00000000137970; P.ptr(0, 2).* = 0.00000010527634;
    P.ptr(1, 0).* = 0.00000000137970; P.ptr(1, 1).* = 0.11247117674074; P.ptr(1, 2).* = 0.00000000702129;
    P.ptr(2, 0).* = 0.00000010527634; P.ptr(2, 1).* = 0.00000000702129; P.ptr(2, 2).* = 0.16745108149516;

    const opt = QuantumDynamicsOptions(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tripleState1D_3"
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

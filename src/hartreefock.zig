//! Module to perform the Hartree-Fock calculation.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const bls = @import("blas.zig"      );
const eig = @import("eigen.zig"     );
const inp = @import("input.zig"     );
const lbt = @import("libint.zig"    );
const lpk = @import("lapack.zig"    );
const mat = @import("matrix.zig"    );
const mth = @import("math.zig"      );
const out = @import("output.zig"    );
const pop = @import("population.zig");
const ten = @import("tensor.zig"    );
const vec = @import("vector.zig"    );

const Basis  = @import("basis.zig" ).Basis ;
const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;
const uncr    = @import("helper.zig").uncr   ;

/// Run the Hartree-Fock calculation with the given options.
pub fn run(comptime T: type, opt: inp.HartreeFockOptions(T), print: bool, allocator: std.mem.Allocator) !out.HartreeFockOutput(T) {
    if (opt.integral.basis == null and (opt.integral.overlap == null or opt.integral.kinetic == null or opt.integral.nuclear == null or opt.integral.coulomb == null)) return error.MissingIntegral;

    if (opt.system_file == null and (opt.system.coords == null or opt.system.atoms == null)) return error.SystemNotFullySpecified;

    if (opt.system.multiplicity != 1) return error.MultiplicityNotImplemented;

    var system: System(T) = undefined; var basis: std.ArrayList(T) = undefined;

    if (opt.system_file == null and opt.system.coords != null) {

        system = System(T){
            .coords = try Matrix(T).init(opt.system.coords.?.len, 3, allocator),
            .atoms  = try Vector(T).init(opt.system.atoms .?.len,    allocator),
            .nocc   = mth.sum(u8, opt.system.atoms.?) / 2,
        };

        for (0..opt.system.atoms.?.len) |i| {
            system.atoms.ptr(i).* = @as(T, @floatFromInt(opt.system.atoms.?[i]));
        }

        for (0..opt.system.coords.?.len) |i| for (0..3) |j| {
            system.coords.ptr(i, j).* = opt.system.coords.?[i][j] * A2AU;
        };
    }

    if (opt.system_file    != null) system = try System(T).read(        opt.system_file.?,    allocator);
    if (opt.integral.basis != null) basis  = try Basis(T).array(system, opt.integral.basis.?, allocator);

    system.nocc -= opt.system.charge / 2; var nbf: usize = 0; var npgs: usize = 0; var mem: f64 = 0; const nocc = system.nocc; const VNN = system.nuclearRepulsion();

    {
        var i: usize = 0; while (opt.integral.basis != null and i < basis.items.len) : (i += 2 * @as(usize, @intFromFloat(basis.items[i])) + 5) {
            const cgs: usize = @as(usize, @intFromFloat((basis.items[i + 1] + 1) * (basis.items[i + 1] + 2))) / 2; nbf += cgs; npgs += @as(usize, @intFromFloat(basis.items[i])) * cgs;
        }

        mem = 8 * @as(f64, @floatFromInt(12 * nbf * nbf + 2 * nbf * nbf * nbf * nbf));
    }

    if (print and opt.integral.basis != null) try std.io.getStdOut().writer().print("\n# OF CONTRACTED GAUSSIAN SHELLS: {d}\n", .{nbf });
    if (print and opt.integral.basis != null) try std.io.getStdOut().writer().print(  "# OF PRIMITIVE  GAUSSIAN SHELLS: {d}\n", .{npgs});

    if (print and opt.integral.basis != null) try std.io.getStdOut().writer().print("\n# OF GB NEEDED TO PERFORM THE HF CALCULATION: {d:.2}\n", .{mem / 1e9});

    var timer = try std.time.Timer.start();

    const S_AO = if (opt.integral.overlap != null) try mat.read(T, opt.integral.overlap.?,    allocator) else try Matrix(T).init(nbf, nbf,                      allocator);
    const T_AO = if (opt.integral.kinetic != null) try mat.read(T, opt.integral.kinetic.?,    allocator) else try Matrix(T).init(nbf, nbf,                      allocator);
    const V_AO = if (opt.integral.nuclear != null) try mat.read(T, opt.integral.nuclear.?,    allocator) else try Matrix(T).init(nbf, nbf,                      allocator);
    const J_AO = if (opt.integral.coulomb != null) try ten.read(T, opt.integral.coulomb.?, 4, allocator) else try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator);

    if (opt.integral.overlap == null) lbt.overlap(S_AO.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));
    if (opt.integral.kinetic == null) lbt.kinetic(T_AO.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));
    if (opt.integral.nuclear == null) lbt.nuclear(V_AO.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));
    if (opt.integral.coulomb == null) lbt.coulomb(J_AO.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));

    nbf = S_AO.rows; var J_AO_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_AO_A.deinit();

    for (0..nbf) |i| for(0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
        J_AO_A.ptr(&[_]usize{i, j, k, l}).* = J_AO.at(&[_]usize{i, j, k, l}) - 0.5 * J_AO.at(&[_]usize{i, l, k, j});
    };

    if (print) try std.io.getStdOut().writer().print("\nINTEGRALS OBTAINED: {}\n", .{std.fmt.fmtDuration(timer.read())});

    var T1 = try Matrix(T).init(nbf, nbf, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nbf, nbf, allocator); defer T2.deinit();
    var T3 = try Matrix(T).init(nbf, nbf, allocator); defer T3.deinit();

    var H_AO = try Matrix(T).init(nbf, nbf, allocator); defer H_AO.deinit();
    var ERR  = try Matrix(T).init(nbf, nbf, allocator); defer  ERR.deinit();

    var DIIS_E = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_E.deinit();
    var DIIS_F = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_F.deinit();

    if (opt.dsize != null) for (0..opt.dsize.?) |_| {
        try DIIS_E.append(try Matrix(T).init(H_AO.rows, H_AO.cols, allocator));
        try DIIS_F.append(try Matrix(T).init(H_AO.rows, H_AO.cols, allocator));
    };

    var F_AO = try Matrix(T).init(nbf, nbf, allocator); F_AO.fill(0);
    var C_MO = try Matrix(T).init(nbf, nbf, allocator); C_MO.fill(0);
    var D_MO = try Matrix(T).init(nbf, nbf, allocator); D_MO.fill(0);
    var E_MO = try Matrix(T).init(nbf, nbf, allocator); E_MO.fill(0);

    var iter: u32 = 0; var E: T = 0; var EP: T = 1; mat.add(T, &H_AO, T_AO, V_AO);

    if (print) try std.io.getStdOut().writer().print("\n{s:4} {s:20} {s:8}\n", .{"ITER", "ENERGY", "|DELTA E|"});

    while (@abs(EP - E) > opt.threshold) : (iter += 1) {

        if (iter == opt.maxiter) return error.MaxIterationsExceeded;

        try eig.contract(&F_AO, D_MO, J_AO_A, &[_]i32{0, 0, 1, 1}); mat.add(T, &F_AO, F_AO, H_AO);

        if (opt.dsize != null and iter > 0) {

            bls.dgemm(&T1, S_AO, false, D_MO, false); bls.dgemm(&T2, T1, false, F_AO, true); bls.dgemm(&T1, F_AO, false, D_MO, false); bls.dgemm(&T3, T1, false, S_AO, true); mat.sub(T, &ERR, T2, T3);

            @memcpy(DIIS_F.items[iter % DIIS_F.items.len].data, F_AO.data);
            @memcpy(DIIS_E.items[iter % DIIS_E.items.len].data,  ERR.data);

            try diisExtrapolate(T, &F_AO, &DIIS_F, &DIIS_E, iter, allocator);
        }

        lpk.dsygvd(&E_MO, &C_MO, F_AO, S_AO, &T1); D_MO.fill(0); EP = E; E = 0;

        for (0..nbf) |i| for (0..nocc) |j| for (0..nbf) |k| {
            D_MO.ptr(i, k).* += 2.0 * C_MO.at(i, j) * C_MO.at(k, j);
        };

        for (0..nbf) |i| for (0..nbf) |j| {
            E += 0.5 * D_MO.at(i, j) * (H_AO.at(i, j) + F_AO.at(i, j));
        };

        if (print) try std.io.getStdOut().writer().print("{d:4} {d:20.14} {e:9.3}\n", .{iter + 1, E + VNN, @abs(EP - E)});
    }

    if (print) try std.io.getStdOut().writer().print("\nHF ENERGY: {d:.14}\n", .{E + VNN});

    const m = try pop.mulliken(T, system, basis, S_AO, D_MO, allocator); defer m.deinit();

    if (print) {
        try std.io.getStdOut().writer().print("\nHF MULLIKEN POPULATION ANALYSIS:\n", .{}); try m.matrix().print(std.io.getStdOut().writer());
    }

    for (0..DIIS_E.items.len) |i| DIIS_E.items[i].deinit();
    for (0..DIIS_F.items.len) |i| DIIS_F.items[i].deinit();

    if (opt.write.coefficient_ao != null) try C_MO.write(opt.write.coefficient_ao.?);
    if (opt.write.coulomb_ao     != null) try J_AO.write(opt.write.coulomb_ao.?    );
    if (opt.write.density_ao     != null) try D_MO.write(opt.write.density_ao.?    );
    if (opt.write.kinetic_ao     != null) try T_AO.write(opt.write.kinetic_ao.?    );
    if (opt.write.nuclear_ao     != null) try V_AO.write(opt.write.nuclear_ao.?    );
    if (opt.write.overlap_ao     != null) try S_AO.write(opt.write.overlap_ao.?    );

    return out.HartreeFockOutput(T){
        .S_AO = S_AO, .T_AO = T_AO, .V_AO = V_AO, .J_AO = J_AO, .C_MO = C_MO, .D_MO = D_MO, .E_MO = E_MO, .F_AO = F_AO, .E = E + VNN, .system = system, .basis = basis
    };
}

/// Extrapolate the DIIS error to obtain a new Fock matrix.
pub fn diisExtrapolate(comptime T: type, F_AO: *Matrix(T), DIIS_F: *std.ArrayList(Matrix(T)), DIIS_E: *std.ArrayList(Matrix(T)), iter: u32, allocator: std.mem.Allocator) !void {
    const size = if (iter < DIIS_F.items.len) iter else DIIS_F.items.len;

    var A   = try Matrix(T  ).init(size + 1, size + 1, allocator); defer   A.deinit();
    var ALU = try Matrix(T  ).init(size + 1, size + 1, allocator); defer ALU.deinit();
    var b   = try Vector(T  ).init(size + 1,           allocator); defer   b.deinit();
    var c   = try Vector(T  ).init(size + 1,           allocator); defer   c.deinit();
    var p   = try Vector(i32).init(size + 1,           allocator); defer   p.deinit();

    A.fill(1); b.fill(0); A.ptr(A.rows - 1, A.cols - 1).* = 0; b.ptr(b.rows - 1).* = 1;

    for (0..size) |i| for (0..size) |j| {

        A.ptr(i, j).* = 0;

        const ii = (iter - size + i + 1) % DIIS_E.items.len;
        const jj = (iter - size + j + 1) % DIIS_E.items.len;

        for (0..DIIS_E.items[0].rows) |k| for (0..DIIS_E.items[0].cols) |l| {
            A.ptr(i, j).* += DIIS_E.items[ii].at(k, l) * DIIS_E.items[jj].at(k, l);
        };
    };

    lpk.dgesv(&c, &ALU, &p, A, b); if (lpk.dgecon(ALU, lpk.dlange(A, '1')) < 1e-12) return else F_AO.fill(0);

    for (0..size) |i| {

        const ii = (iter - size + i + 1) % DIIS_E.items.len;

        for (0..F_AO.rows) |j| for (0..F_AO.cols) |k| {
            F_AO.ptr(j, k).* += c.at(i) * DIIS_F.items[ii].at(j, k);
        };
    }
}

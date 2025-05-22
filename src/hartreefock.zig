//! Module to perform the Hartree-Fock calculation.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const bls = @import("blas.zig"      );
const edf = @import("energydiff.zig");
const eig = @import("eigen.zig"     );
const inp = @import("input.zig"     );
const lbt = @import("libint.zig"    );
const lpk = @import("lapack.zig"    );
const mat = @import("matrix.zig"    );
const mth = @import("math.zig"      );
const opm = @import("optimize.zig"  );
const out = @import("output.zig"    );
const pop = @import("population.zig");
const prp = @import("property.zig"  );
const sys = @import("system.zig"    );
const ten = @import("tensor.zig"    );
const tns = @import("transform.zig" );
const vec = @import("vector.zig"    );

const Basis  = @import("basis.zig" ).Basis ;
const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;
const uncr    = @import("helper.zig").uncr   ;

/// Main function to run the Hartree-Fock calculation with the given options.
pub fn run(comptime T: type, opt: inp.HartreeFockOptions(T), print: bool, allocator: std.mem.Allocator) !out.HartreeFockOutput(T) {
     var system = try sys.load(T, opt.system, opt.system_file, allocator); defer system.deinit();

    if (print) {
        try std.io.getStdOut().writer().print("\nSYSTEM:\n", .{}); try system.print(std.io.getStdOut().writer());
    }

    return runFull(T, opt, &system, print, allocator);
}

/// Secondary function to run the Hartree-Fock calculation with the given options from a reference geometry.
pub fn runFull(comptime T: type, opt: inp.HartreeFockOptions(T), system: *System(T), print: bool, allocator: std.mem.Allocator) !out.HartreeFockOutput(T) {
    try checkErrors(T, opt); const name = try std.fmt.allocPrint(allocator, "HF", .{}); defer allocator.free(name);

    if (opt.optimize != null) {
        const optsys = try opm.optimize(T, opt, system.*, hfFull, name, print, allocator); system.deinit(); system.* = optsys;
    }

    if (print) {
        if (opt.optimize != null) {try std.io.getStdOut().writer().print("\n{s} OPTIMIZED SYSTEM:\n", .{name}); try system.print(std.io.getStdOut().writer());}
    }

    var hf = try hfFull(T, opt, system.*, print, allocator);

    if (opt.mulliken) hf.mulliken = try pop.mulliken(T, system.*, hf.basis, hf.S_AS, hf.D_AS, allocator);

    if (print) {
        if (hf.mulliken != null) {try std.io.getStdOut().writer().print("\n{s} MULLIKEN POPULATION ANALYSIS:\n", .{name}); try hf.mulliken.?.matrix().print(std.io.getStdOut().writer());}
    }

    if (opt.gradient != null) hf.G = try edf.gradient(T, opt, system.*, hfFull, name, true, allocator);

    if (print) {
        if (hf.G != null) {try std.io.getStdOut().writer().print("\n{s} GRADIENT:\n", .{name}); try hf.G.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null) hf.H = try edf.hessian(T, opt, system.*, hfFull, name, true, allocator);

    if (print and opt.hessian != null and opt.hessian.?.print) {
        if (hf.H != null) {try std.io.getStdOut().writer().print("\n{s} HESSIAN:\n", .{name}); try hf.H.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null and opt.hessian.?.freq) hf.freqs = try prp.freq(T, system.*, hf.H.?, allocator);

    if (print) {
        if (hf.freqs != null) {try std.io.getStdOut().writer().print("\n{s} HARMONIC FREQUENCIES:\n", .{name}); try hf.freqs.?.matrix().print(std.io.getStdOut().writer());}
    }

    return hf;
}

/// Function to actually run the Hartree-Fock energy calculation on the provided system without additional calculation.
pub fn hfFull(comptime T: type, opt: inp.HartreeFockOptions(T), system: System(T), print: bool, allocator: std.mem.Allocator) !out.HartreeFockOutput(T) {
    var basis: std.ArrayList(T) = undefined; if (opt.integral.basis != null) basis = try Basis(T).array(system, opt.integral.basis.?, allocator);

    var nbf: usize = 0; var npgs: usize = 0; var mem: f64 = 0; const VNN = system.nuclearRepulsion();

    {
        var i: usize = 0; while (opt.integral.basis != null and i < basis.items.len) : (i += 2 * @as(usize, @intFromFloat(basis.items[i])) + 5) {
            const cgs: usize = @as(usize, @intFromFloat((basis.items[i + 1] + 1) * (basis.items[i + 1] + 2))) / 2; nbf += cgs; npgs += @as(usize, @intFromFloat(basis.items[i])) * cgs;
        }

        nbf = if (opt.generalized) 2 * nbf else nbf; npgs = if (opt.generalized) 2 * npgs else npgs;

        mem = 8 * @as(f64, @floatFromInt(12 * nbf * nbf + 2 * nbf * nbf * nbf * nbf));
    }

    if (print) try std.io.getStdOut().writer().print("\n# NUMBER OF ELECTRONS: {d}\n", .{system.getElectrons()});

    if (print and opt.integral.basis != null) try std.io.getStdOut().writer().print("\n# OF CONTRACTED GAUSSIAN SHELLS: {d}\n", .{nbf });
    if (print and opt.integral.basis != null) try std.io.getStdOut().writer().print(  "# OF PRIMITIVE  GAUSSIAN SHELLS: {d}\n", .{npgs});

    if (print and opt.integral.basis != null) try std.io.getStdOut().writer().print("\n# OF GB NEEDED TO PERFORM THE HF CALCULATION: {d:.2}\n", .{mem / 1e9});

    var timer = try std.time.Timer.start(); const nbf_spatial = if (opt.generalized) nbf / 2 else nbf;

    var S_A = if (opt.integral.overlap != null) try mat.read(T, opt.integral.overlap.?,    allocator) else try Matrix(T).init(nbf_spatial, nbf_spatial,                                      allocator);
    var T_A = if (opt.integral.kinetic != null) try mat.read(T, opt.integral.kinetic.?,    allocator) else try Matrix(T).init(nbf_spatial, nbf_spatial,                                      allocator);
    var V_A = if (opt.integral.nuclear != null) try mat.read(T, opt.integral.nuclear.?,    allocator) else try Matrix(T).init(nbf_spatial, nbf_spatial,                                      allocator);
    var J_A = if (opt.integral.coulomb != null) try ten.read(T, opt.integral.coulomb.?, 4, allocator) else try Tensor(T).init(&[_]usize{nbf_spatial, nbf_spatial, nbf_spatial, nbf_spatial}, allocator);

    if (opt.integral.overlap == null) lbt.overlap(S_A.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));
    if (opt.integral.kinetic == null) lbt.kinetic(T_A.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));
    if (opt.integral.nuclear == null) lbt.nuclear(V_A.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));
    if (opt.integral.coulomb == null) lbt.coulomb(J_A.data, system, try Basis(f64).array(system, opt.integral.basis.?, std.heap.page_allocator));

    if (opt.generalized) {

        var S_AS = try Matrix(T).init(nbf, nbf,                      allocator);
        var T_AS = try Matrix(T).init(nbf, nbf,                      allocator);
        var V_AS = try Matrix(T).init(nbf, nbf,                      allocator);
        var J_AS = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator);

        tns.oneAO2AS(T, &S_AS, S_A); tns.oneAO2AS(T, &T_AS, T_A); tns.oneAO2AS(T, &V_AS, V_A); tns.twoAO2AS(T, &J_AS, J_A);

        S_A.deinit(); S_A = S_AS;
        T_A.deinit(); T_A = T_AS;
        V_A.deinit(); V_A = V_AS;
        J_A.deinit(); J_A = J_AS;
    }

    const nocc: usize = if (opt.generalized) system.getElectrons() else system.getElectrons() / 2;

    var J_A_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator);

    const jf: T = if (opt.generalized) 1 else 0.5; for (0..nbf) |i| for(0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
        J_A_A.ptr(&[_]usize{i, j, k, l}).* = J_A.at(&[_]usize{i, j, k, l}) - jf * J_A.at(&[_]usize{i, l, k, j});
    };

    if (print) try std.io.getStdOut().writer().print("\nINTEGRALS OBTAINED: {}\n", .{std.fmt.fmtDuration(timer.read())});

    var T1 = try Matrix(T).init(nbf, nbf, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nbf, nbf, allocator); defer T2.deinit();
    var T3 = try Matrix(T).init(nbf, nbf, allocator); defer T3.deinit();

    var H_A = try Matrix(T).init(nbf, nbf, allocator); defer H_A.deinit();
    var ERR = try Matrix(T).init(nbf, nbf, allocator); defer ERR.deinit();

    var DIIS_E = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_E.deinit();
    var DIIS_F = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_F.deinit();

    if (opt.dsize != null) for (0..opt.dsize.?) |_| {
        try DIIS_E.append(try Matrix(T).init(H_A.rows, H_A.cols, allocator));
        try DIIS_F.append(try Matrix(T).init(H_A.rows, H_A.cols, allocator));
    };

    var F_A = try Matrix(T).init(nbf, nbf, allocator); F_A.fill(0);
    var C_A = try Matrix(T).init(nbf, nbf, allocator); C_A.fill(0);
    var D_A = try Matrix(T).init(nbf, nbf, allocator); D_A.fill(0);
    var E_M = try Matrix(T).init(nbf, nbf, allocator); E_M.fill(0);

    var iter: u32 = 0; var E: T = 0; var EP: T = 1; mat.add(T, &H_A, T_A, V_A);

    if (print) try std.io.getStdOut().writer().print("\nHF SELF-CONSISTENT FIELD:\n{s:4} {s:20} {s:8}\n", .{"ITER", "ENERGY", "|DELTA E|"});

    while (@abs(EP - E) > opt.threshold) : (iter += 1) {

        if (iter == opt.maxiter) return error.MaxIterationsExceeded;

        try eig.contract(&F_A, D_A, J_A_A, &[_]i32{0, 0, 1, 1}); mat.add(T, &F_A, F_A, H_A);

        if (opt.dsize != null and iter > 0) {

            bls.dgemm(&T1, S_A, false, D_A, false); bls.dgemm(&T2, T1, false, F_A, true); bls.dgemm(&T1, F_A, false, D_A, false); bls.dgemm(&T3, T1, false, S_A, true); mat.sub(T, &ERR, T2, T3);

            @memcpy(DIIS_F.items[iter % DIIS_F.items.len].data, F_A.data);
            @memcpy(DIIS_E.items[iter % DIIS_E.items.len].data,  ERR.data);

            try diisExtrapolate(T, &F_A, &DIIS_F, &DIIS_E, iter, allocator);
        }

        lpk.dsygvd(&E_M, &C_A, F_A, S_A, &T1); D_A.fill(0); EP = E; E = 0;

        const df: T = if(opt.generalized) 1 else 2; for (0..nbf) |i| for (0..nocc) |j| for (0..nbf) |k| {
            D_A.ptr(i, k).* += df * C_A.at(i, j) * C_A.at(k, j);
        };

        for (0..nbf) |i| for (0..nbf) |j| {
            E += 0.5 * D_A.at(i, j) * (H_A.at(i, j) + F_A.at(i, j));
        };

        if (print) try std.io.getStdOut().writer().print("{d:4} {d:20.14} {e:9.3}\n", .{iter + 1, E + VNN, @abs(EP - E)});
    }

    if (print) try std.io.getStdOut().writer().print("\nHF ENERGY: {d:.14}\n", .{E + VNN});

    for (0..DIIS_E.items.len) |i| DIIS_E.items[i].deinit();
    for (0..DIIS_F.items.len) |i| DIIS_F.items[i].deinit();

    J_A_A.deinit();

    if (opt.write.coefficient_ao != null) try C_A.write(opt.write.coefficient_ao.?);
    if (opt.write.coulomb_ao     != null) try J_A.write(opt.write.coulomb_ao.?    );
    if (opt.write.density_ao     != null) try D_A.write(opt.write.density_ao.?    );
    if (opt.write.kinetic_ao     != null) try T_A.write(opt.write.kinetic_ao.?    );
    if (opt.write.nuclear_ao     != null) try V_A.write(opt.write.nuclear_ao.?    );
    if (opt.write.overlap_ao     != null) try S_A.write(opt.write.overlap_ao.?    );

    if (!opt.generalized) {

        var S_AS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var T_AS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var V_AS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var F_AS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var C_AS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var D_AS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var E_MS = try Matrix(T).init(2 * nbf, 2 * nbf,                              allocator);
        var J_AS = try Tensor(T).init(&[_]usize{2 * nbf, 2 * nbf, 2 * nbf, 2 * nbf}, allocator);

        tns.oneAO2AS(T, &S_AS, S_A); tns.oneAO2AS(T, &T_AS, T_A); tns.oneAO2AS(T, &V_AS, V_A); tns.oneAO2AS(T, &F_AS, F_A); tns.oneAO2AS(T, &D_AS, D_A);

        tns.cfsAO2AS(T, &C_AS, C_A); tns.twoAO2AS(T, &J_AS, J_A); E_MS.fill(0); for (0..2 * nbf) |i| E_MS.ptr(i, i).* = E_M.at(i / 2, i / 2);

        S_A.deinit(); S_A = S_AS;
        T_A.deinit(); T_A = T_AS;
        V_A.deinit(); V_A = V_AS;
        F_A.deinit(); F_A = F_AS;
        C_A.deinit(); C_A = C_AS;
        D_A.deinit(); D_A = D_AS;
        E_M.deinit(); E_M = E_MS;
        J_A.deinit(); J_A = J_AS;

        mat.muls(T, &D_A, D_A, 0.5);
    }

    return out.HartreeFockOutput(T){
        .S_AS = S_A, .T_AS = T_A, .V_AS = V_A, .J_AS = J_A, .C_AS = C_A, .D_AS = D_A, .E_MS = E_M, .F_AS = F_A, .E = E + VNN, .mulliken = null, .G = null, .H = null, .freqs = null, .basis = basis
    };
}

/// Check the errors int the provided input.
pub fn checkErrors(comptime T: type, opt: inp.HartreeFockOptions(T)) !void {
    const calcInt = opt.integral.overlap == null or opt.integral.kinetic == null or opt.integral.nuclear == null or opt.integral.coulomb == null;
    const readInt = opt.integral.overlap != null or opt.integral.kinetic != null or opt.integral.nuclear != null or opt.integral.coulomb != null;

    if (opt.integral.basis == null and calcInt) return error.MissingIntegral;

    if (opt.gradient != null and readInt) return error.CantUseGradientWithProvidedIntegrals;
    if (opt.hessian  != null and readInt) return  error.CantUseHessianWithProvidedIntegrals;

    if (!opt.generalized and @rem(@abs(opt.system.charge), 2) == 1) return error.UseGeneralizedHartreeFockForOddCharge;

    if (opt.dsize != null and opt.dsize.? < 1) return error.InvalidDIISSize;
}

/// Extrapolate the DIIS error to obtain a new Fock matrix.
pub fn diisExtrapolate(comptime T: type, F_A: *Matrix(T), DIIS_F: *std.ArrayList(Matrix(T)), DIIS_E: *std.ArrayList(Matrix(T)), iter: u32, allocator: std.mem.Allocator) !void {
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

    lpk.dgesv(&c, &ALU, &p, A, b); if (lpk.dgecon(ALU, lpk.dlange(A, '1')) < 1e-12) return else F_A.fill(0);

    for (0..size) |i| {

        const ii = (iter - size + i + 1) % DIIS_E.items.len;

        for (0..F_A.rows) |j| for (0..F_A.cols) |k| {
            F_A.ptr(j, k).* += c.at(i) * DIIS_F.items[ii].at(j, k);
        };
    }
}

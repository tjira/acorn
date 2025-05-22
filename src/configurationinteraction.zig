//! Module for Configuration Interaction (CI) calculations.

const std = @import("std");

const cas = @import("completeactivespace.zig");
const edf = @import("energydiff.zig"         );
const hfm = @import("hartreefock.zig"        );
const inp = @import("input.zig"              );
const lpk = @import("lapack.zig"             );
const mat = @import("matrix.zig"             );
const mth = @import("math.zig"               );
const opm = @import("optimize.zig"           );
const out = @import("output.zig"             );
const prp = @import("property.zig"           );
const sys = @import("system.zig"             );
const ten = @import("tensor.zig"             );
const tns = @import("transform.zig"          );

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

/// Main function to run the CI calculation with the given options.
pub fn run(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    var system = try sys.load(T, opt.hartree_fock.system, opt.hartree_fock.system_file, allocator); defer system.deinit();

    if (print) {
        try std.io.getStdOut().writer().print("\nSYSTEM:\n", .{}); try system.print(std.io.getStdOut().writer());
    }

    return runFull(T, opt, &system, print, allocator);
}

/// Secondary function to run the CI calculation with the given options from a reference geometry.
pub fn runFull(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), system: *System(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    try checkErrors(T, opt); const name = try std.fmt.allocPrint(allocator, "CI", .{}); defer allocator.free(name);

    if (opt.optimize != null) {
        const optsystem = try opm.optimize(T, opt, system.*, ciFull, name, print, allocator); system.deinit(); system.* = optsystem;
    }

    if (print) {
        if (opt.optimize != null) {try std.io.getStdOut().writer().print("\n{s} OPTIMIZED SYSTEM:\n", .{name}); try system.print(std.io.getStdOut().writer());}
    }

    const hf = try hfm.runFull(T, opt.hartree_fock, system, print, allocator); var ci = try ciPost(T, opt, system.*, hf, print, allocator);

    if (opt.gradient != null) ci.G = try edf.gradient(T, opt, system.*, ciFull, name, true, allocator);

    if (print) {
        if (ci.G != null) {try std.io.getStdOut().writer().print("\n{s} GRADIENT:\n", .{name}); try ci.G.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null) ci.H = try edf.hessian(T, opt, system.*, ciFull, name, true, allocator);

    if (print and opt.hessian != null and opt.hessian.?.print) {
        if (ci.H != null) {try std.io.getStdOut().writer().print("\n{s} HESSIAN:\n", .{name}); try ci.H.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null and opt.hessian.?.freq) ci.freqs = try prp.freq(T, system.*, ci.H.?, allocator);

    if (print) {
        if (ci.freqs != null) {try std.io.getStdOut().writer().print("\n{s} HARMONIC FREQUENCIES:\n", .{name}); try ci.freqs.?.matrix().print(std.io.getStdOut().writer());}
    }

    return ci;
}

/// Function to actually run the CI energy calculation on the provided system without additional calculation.
pub fn ciFull(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), system: System(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    const hf = try hfm.hfFull(T, opt.hartree_fock, system, print, allocator); return ciPost(T, opt, system, hf, print, allocator);
}

/// Function to run the CI energy calculation on the provided system with Hartree-Fock output.
pub fn ciPost(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), system: System(T), hf: out.HartreeFockOutput(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    const nbf = hf.S_AS.rows; const nocc = system.getElectrons(); var AS = [2]usize{nocc, nbf};

    if (opt.active_space != null) {
        AS[0] = opt.active_space.?[0];
        AS[1] = opt.active_space.?[1];
    }

    if (AS[0] > nocc or AS[1] > nbf or AS[0] > AS[1] or AS[1] - AS[0] > nbf - nocc) return error.InvalidActiveSpace;

    if (print) try std.io.getStdOut().writer().print("\nCI ACTIVE SPACE: {d} ELECTRONS IN {d} SPINORBITALS\n", .{AS[0], AS[1]});

    var D = try cas.generateCasDeterminants(AS[0], AS[1], nocc, opt.hartree_fock.generalized, allocator); defer D.deinit();

    var H_MS   = try Matrix(T).init(nbf, nbf,                      allocator); defer   H_MS.deinit();
    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();

    try transform(T, &H_MS, &J_MS_A, hf.T_AS, hf.V_AS, hf.J_AS, hf.C_AS, allocator);

    if (print) try std.io.getStdOut().writer().print("\nNUMBER OF CI DETERMINANTS: {d}\n", .{D.rows});

    var A = try Vector(usize).init(nocc,           allocator); defer A.deinit();
    var H = try Matrix(T    ).init(D.rows, D.rows, allocator); defer H.deinit();
    var E = try Matrix(T    ).init(D.rows, D.rows, allocator); defer E.deinit();
    var C = try Matrix(T    ).init(D.rows, D.rows, allocator); defer C.deinit();

    H.fill(0);

    for (0..D.rows) |i| for (i..D.rows) |j| {

        var so: [4]usize = undefined; var diff: u32 = 0; var k: usize = 0;

        const sign = try alignDeterminant(&A, D.row(i).vector(), D.row(j).vector());

        for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {diff += 1;};

        if (diff > 2) continue;

        if (diff == 1) for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {
            so[0] = A.at(l); so[1] = D.at(i, l);
        };

        if (diff == 2) for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {
            so[if (k == 0) 0 else 1] = A.at(l); so[if (k == 0) 2 else 3] = D.at(i, l); k += 1;
        };

        H.ptr(i, j).* = asfloat(T, sign) * try slater(T, A, so[0..diff * 2], H_MS, J_MS_A); H.ptr(j, i).* = H.at(i, j);
    };

    lpk.dsyevd(&E, &C, H);

    if (print) try std.io.getStdOut().writer().print("\nCI ENERGY: {d:.14}\n", .{E.at(0, 0) + system.nuclearRepulsion()});

    return .{
        .hf = hf, .E = E.at(0, 0) + system.nuclearRepulsion(), .G = null, .H = null, .freqs = null,
    };
}

/// Check for errors in the CI options.
pub fn checkErrors(comptime T: type, opt: inp.ConfigurationInteractionOptions(T)) !void {
    const readInt = opt.hartree_fock.integral.overlap != null or opt.hartree_fock.integral.kinetic != null or opt.hartree_fock.integral.nuclear != null or opt.hartree_fock.integral.coulomb != null;

    if (opt.gradient != null and readInt) return error.CantUseGradientWithProvidedIntegrals;
    if (opt.hessian  != null and readInt) return  error.CantUseHessianWithProvidedIntegrals;

    if (!opt.hartree_fock.generalized and opt.active_space != null and opt.active_space.?[1] % 2 == 1) return error.UseGeneralizedHartreeFockForOddActiveOrbitals;

    if (opt.hartree_fock.optimize != null) return error.CantOptimizeOnHartreeFockAndUseConfigurationInteraction;
}

/// Aligns the vector C to the vector B. The result is stored in the vector A and the sign of the permutation is returned.
fn alignDeterminant(A: *Vector(usize), B: Vector(usize), C: Vector(usize)) !i32 {
    @memcpy(A.data, C.data); var k: i32 = 0; var sign: i32 = 1;

    while (k < A.rows) : (k += 1) {
        if (A.at(@intCast(k)) != B.at(@intCast(k))) for (@as(usize, @intCast(k)) + 1..A.rows) |l| if (A.at(@intCast(k)) == B.at(l) or A.at(l) == B.at(@intCast(k))) {
            std.mem.swap(usize, A.ptr(@intCast(k)), A.ptr(l)); sign *= -1; k -= 1; break;
        };
    }

    return sign;
}

/// Slater-Condon rules for the CI calculations.
fn slater(comptime T: type, A: Vector(usize), so: []const usize, H_MS: Matrix(T), J_MS_A: Tensor(T)) !T {
    var hij: T = 0;

    if (so.len / 2 == 0) {

        for (0..A.rows) |l| {
            hij += H_MS.at(A.at(l), A.at(l));
        }

        for (0..A.rows) |l| for (0..A.rows) |m| {
            hij += 0.5 * J_MS_A.at(&[_]usize{A.at(l), A.at(m), A.at(l), A.at(m)});
        };
    }

    else if (so.len / 2 == 1) {

        hij += H_MS.at(so[0], so[1]);

        for (0..A.rows) |m| if (A.at(m) != so[0]) {
            hij += J_MS_A.at(&[_]usize{so[0], A.at(m), so[1], A.at(m)});
        };
    }

    else if (so.len / 2 == 2) {
        hij = J_MS_A.at(&[_]usize{so[0], so[1], so[2], so[3]});
    }
    
    return hij;
}

/// Function to perform all integrals transformations used in the CI calculations.
fn transform(comptime T: type, H_MS: *Matrix(T), J_MS_A: *Tensor(T), T_AS: Matrix(T), V_AS: Matrix(T), J_AS: Tensor(T), C_AS: Matrix(T), allocator: std.mem.Allocator) !void {
    var H_AS = try Matrix(T).init(T_AS.rows, T_AS.cols,                                                          allocator); defer H_AS.deinit();
    var J_MS = try Tensor(T).init(&[_]usize{J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0]}, allocator); defer J_MS.deinit();

    mat.add(T, &H_AS, T_AS, V_AS); tns.oneAO2MO(T, H_MS, H_AS, C_AS); try tns.twoAO2MO(T, &J_MS, J_AS, C_AS);

    for (0..J_MS.shape[0]) |i| for (0..J_MS.shape[1]) |j| for (0..J_MS.shape[2]) |k| for (0..J_MS.shape[3]) |l| {
        J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS.at(&[_]usize{i, j, k, l}) - J_MS.at(&[_]usize{i, l, k, j});
    };
}

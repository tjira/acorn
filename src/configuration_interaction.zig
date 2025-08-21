//! Module for Configuration Interaction (CI) calculations.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const cas = @import("complete_active_space.zig");
const edf = @import("energy_derivative.zig"    );
const hfm = @import("hartree_fock.zig"         );
const hlp = @import("helper.zig"               );
const inp = @import("input.zig"                );
const mat = @import("matrix.zig"               );
const mth = @import("math.zig"                 );
const opm = @import("optimize.zig"             );
const out = @import("output.zig"               );
const prp = @import("property.zig"             );
const sys = @import("system.zig"               );
const ten = @import("tensor.zig"               );
const tns = @import("transform.zig"            );

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Main function to run the CI calculation with the given options.
pub fn run(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    var system = try sys.load(T, opt.hartree_fock.system, opt.hartree_fock.system_file, allocator); defer system.deinit();

    if (print) {
        try hlp.print(std.fs.File.stdout(), "\nSYSTEM:\n", .{}); try system.print(std.fs.File.stdout());
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
        if (opt.optimize != null) {try hlp.print(std.fs.File.stdout(), "\n{s} OPTIMIZED SYSTEM:\n", .{name}); try system.print(std.fs.File.stdout());}
    }

    const hf = try hfm.runFull(T, opt.hartree_fock, system, print, allocator); var ci = try ciPost(T, opt, system.*, hf, print, allocator);

    if (opt.gradient != null) ci.G = try edf.gradient(T, opt, system.*, ciFull, name, true, allocator);

    if (print) {
        if (ci.G != null) {try hlp.print(std.fs.File.stdout(), "\n{s} GRADIENT:\n", .{name}); try ci.G.?.print(std.fs.File.stdout());}
    }

    if (opt.hessian != null) ci.H = try edf.hessian(T, opt, system.*, ciFull, name, true, allocator);

    if (print and opt.hessian != null and opt.print.hessian) {
        if (ci.H != null) {try hlp.print(std.fs.File.stdout(), "\n{s} HESSIAN:\n", .{name}); try ci.H.?.print(std.fs.File.stdout());}
    }

    if (opt.hessian != null and opt.hessian.?.freq) ci.freqs = try prp.freq(T, system.*, ci.H.?, allocator);

    if (print) {
        if (ci.freqs != null) {try hlp.print(std.fs.File.stdout(), "\n{s} HARMONIC FREQUENCIES:\n", .{name}); try ci.freqs.?.matrix().print(std.fs.File.stdout());}
    }

    if (opt.gradient != null and opt.write.gradient != null) try ci.G.?.write(opt.write.gradient.?);
    if (opt.hessian  != null and opt.write.hessian  != null) try ci.H.?.write(opt.write.hessian.? );

    return ci;
}

/// Function to actually run the CI energy calculation on the provided system without additional calculation.
pub fn ciFull(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), system: System(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    const hf = try hfm.hfFull(T, opt.hartree_fock, system, print, allocator); return ciPost(T, opt, system, hf, print, allocator);
}

/// Function to run the CI energy calculation on the provided system with Hartree-Fock output.
pub fn ciPost(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), system: System(T), hf: out.HartreeFockOutput(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    const nbf = if (opt.hartree_fock.generalized) hf.S_A.rows else 2 * hf.S_A.rows; const nocc = system.getElectrons(); var AS = [2]usize{nocc, nbf};

    if (opt.active_space != null) {
        AS[0] = opt.active_space.?[0];
        AS[1] = opt.active_space.?[1];
    }

    if (AS[0] > nocc or AS[1] > nbf or AS[0] > AS[1] or AS[1] - AS[0] > nbf - nocc) return error.InvalidActiveSpace;

    if (print) try hlp.print(std.fs.File.stdout(), "\nCI ACTIVE SPACE: {d} ELECTRONS IN {d} SPINORBITALS\n", .{AS[0], AS[1]});

    var D = try cas.generateCasDeterminants(AS[0], AS[1], nocc, opt.hartree_fock.generalized, allocator); defer D.deinit();

    var H_MS   = try Matrix(T).init(nbf, nbf,                      allocator); defer   H_MS.deinit();
    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();

    try transform(T, &H_MS, &J_MS_A, hf.T_A, hf.V_A, hf.J_A, hf.C_A, allocator);

    if (print) try hlp.print(std.fs.File.stdout(), "\nNUMBER OF CI DETERMINANTS: {d}\n", .{D.rows});

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

        H.ptr(i, j).* = hlp.asfloat(T, sign) * try slater(T, A, so[0..diff * 2], H_MS, J_MS_A); H.ptr(j, i).* = H.at(i, j);
    };

    try cwp.Lapack(T).dsyevd(&E, &C, H);

    if (print) try hlp.print(std.fs.File.stdout(), "\nCI ENERGY: {d:.14}\n", .{E.at(0, 0) + system.nuclearRepulsion()});

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
    C.memcpy(A.*); var k: i32 = 0; var sign: i32 = 1;

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
fn transform(comptime T: type, H_MS: *Matrix(T), J_MS_A: *Tensor(T), T_A: Matrix(T), V_A: Matrix(T), J_A: ?Tensor(T), C_A: Matrix(T), allocator: std.mem.Allocator) !void {
    if (J_A == null) return error.ConfigurationInteractionRequiresCoulombIntegrals;

    var H_A = try Matrix(T).init(T_A.rows, T_A.cols, allocator); defer H_A.deinit(); mat.add(T, &H_A, T_A, V_A);

    if (T_A.rows != H_MS.rows) {try tns.oneAO2MS(T, H_MS,   H_A,   C_A           );} else {try tns.oneAO2MO(T, H_MS,   H_A,   C_A           );}
    if (T_A.rows != H_MS.rows) {try tns.twoAO2MS(T, J_MS_A, J_A.?, C_A, allocator);} else {try tns.twoAO2MO(T, J_MS_A, J_A.?, C_A, allocator);}

    for (0..J_MS_A.shape[0]) |i| for (0..J_MS_A.shape[1]) |j| for (0..J_MS_A.shape[2]) |k| for (0..J_MS_A.shape[3]) |l| {
        J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS_A.at(&[_]usize{i, j, k, l}) - J_MS_A.at(&[_]usize{i, l, k, j});
    };
}

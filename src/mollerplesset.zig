//! Module to perform the Moller-Plesset calculations.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const edf = @import("energydiff.zig" );
const hfm = @import("hartreefock.zig");
const inp = @import("input.zig"      );
const out = @import("output.zig"     );
const prp = @import("property.zig"   );
const sys = @import("system.zig"     );
const ten = @import("tensor.zig"     );
const tns = @import("transform.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Main function to run the Moller-Plesset calculations.
pub fn run(comptime T: type, opt: inp.MollerPlessetOptions(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    try checkErrors(T, opt); const system = try sys.load(T, opt.hartree_fock.system, opt.hartree_fock.system_file, allocator); defer system.deinit();

    const hf = try hfm.run(T, opt.hartree_fock, print, allocator); var mp = try mpPost(T, opt, system, hf, print, allocator);

    if (opt.gradient != null) mp.G = try edf.gradient(T, opt, system, mpFull, allocator);

    if (print) {
        if (mp.G != null) {try std.io.getStdOut().writer().print("\nMP{d} GRADIENT:\n", .{opt.order}); try mp.G.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null) mp.H = try edf.hessian(T, opt, system, mpFull, allocator);

    if (print) {
        if (mp.H != null) {try std.io.getStdOut().writer().print("\nMP{d} HESSIAN:\n", .{opt.order}); try mp.H.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null and opt.hessian.?.freq) mp.f = try prp.freq(T, system, mp.H.?, allocator);

    if (print) {
        if (mp.f != null) {try std.io.getStdOut().writer().print("\nMP{d} HARMONIC FREQUENCIES:\n", .{opt.order}); try mp.f.?.matrix().print(std.io.getStdOut().writer());}
    }

    return mp;
}

/// Performs the MP calculation with the given options.
pub fn mpFull(comptime T: type, opt: inp.MollerPlessetOptions(T), system: System(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    const hf = try hfm.hfFull(T, opt.hartree_fock, system, print, allocator); return mpPost(T, opt, system, hf, print, allocator);
}

/// Performs the MP calculation with the given options, provided the Hartree-Fock calculation is already done.
pub fn mpPost(comptime T: type, opt: inp.MollerPlessetOptions(T), system: System(T), hf: out.HartreeFockOutput(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    const nbf = hf.S_AS.rows; const nocc = system.getElectrons();

    var F_MS   = try Matrix(T).init(nbf, nbf,                      allocator); defer   F_MS.deinit();
    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();

    try transform(T, &F_MS, &J_MS_A, hf.F_AS, hf.J_AS, hf.C_AS, allocator);

    var E: T = 0;

    if (opt.order >= 2) E += mp2(T, F_MS, J_MS_A, nocc);

    if (print) try std.io.getStdOut().writer().print("\nMP{d} ENERGY: {d:.14}\n", .{opt.order, hf.E + E});

    return .{
        .hf = hf, .E = hf.E + E, .G = null, .H = null, .f = null,
    };
}

/// Check for errors in the Moller-Plesset options.
pub fn checkErrors(comptime T: type, opt: inp.MollerPlessetOptions(T)) !void {
    const readInt = opt.hartree_fock.integral.overlap != null or opt.hartree_fock.integral.kinetic != null or opt.hartree_fock.integral.nuclear != null or opt.hartree_fock.integral.coulomb != null;

    if (opt.gradient != null and readInt) return error.CantUseGradientWithProvidedIntegrals;
    if (opt.hessian  != null and readInt) return  error.CantUseHessianWithProvidedIntegrals;

    if (opt.order != 2) return error.PerturbationOrderNotImplemented;
}

/// Returns the second-order Moller-Plesset energy.
pub fn mp2(comptime T: type, F_MS: Matrix(T), J_MS_A: Tensor(T), nocc: usize) T {
    var E: T = 0;

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..J_MS_A.shape[0]) |a| for (nocc..J_MS_A.shape[0]) |b| {
        E += 0.25 * J_MS_A.at(&[_]usize{i, j, a, b}) * J_MS_A.at(&[_]usize{a, b, i, j}) / (F_MS.at(i, i) + F_MS.at(j, j) - F_MS.at(a, a) - F_MS.at(b, b));
    };

    return E;
}

/// Function to perform all integrals transformations used in the Moller-Plesset calculations.
pub fn transform(comptime T: type, F_MS: *Matrix(T), J_MS_A: *Tensor(T), F_AS: Matrix(T), J_AS: Tensor(T), C_AS: Matrix(T), allocator: std.mem.Allocator) !void {
    var J_MS = try Tensor(T).init(&[_]usize{J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0]}, allocator); defer J_MS.deinit();

    tns.oneAO2MO(T, F_MS, F_AS, C_AS); try tns.twoAO2MO(T, &J_MS, J_AS, C_AS);

    for (0..J_MS.shape[0]) |i| for (0..J_MS.shape[1]) |j| for (0..J_MS.shape[2]) |k| for (0..J_MS.shape[3]) |l| {
        J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS.at(&[_]usize{i, j, k, l}) - J_MS.at(&[_]usize{i, l, k, j});
    };
}

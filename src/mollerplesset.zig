//! Module to perform the Moller-Plesset calculations.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const edf = @import("energydiff.zig" );
const hfm = @import("hartreefock.zig");
const inp = @import("input.zig"      );
const opm = @import("optimize.zig"   );
const out = @import("output.zig"     );
const prp = @import("property.zig"   );
const sys = @import("system.zig"     );
const ten = @import("tensor.zig"     );
const tns = @import("transform.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Main function to run the Moller-Plesset calculation with the given options.
pub fn run(comptime T: type, opt: inp.MollerPlessetOptions(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    var system = try sys.load(T, opt.hartree_fock.system, opt.hartree_fock.system_file, allocator); defer system.deinit();

    if (print) {
        try std.io.getStdOut().writer().print("\nSYSTEM:\n", .{}); try system.print(std.io.getStdOut().writer());
    }

    return runFull(T, opt, &system, print, allocator);
}

/// Secondary function to run the Moller-Plesset calculation with the given options from a reference geometry.
pub fn runFull(comptime T: type, opt: inp.MollerPlessetOptions(T), system: *System(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    try checkErrors(T, opt); const name = try std.fmt.allocPrint(allocator, "MP{d}", .{opt.order}); defer allocator.free(name);

    if (opt.optimize != null) {
        const optsystem = try opm.optimize(T, opt, system.*, mpFull, name, print, allocator); system.deinit(); system.* = optsystem;
    }

    if (print) {
        if (opt.optimize != null) {try std.io.getStdOut().writer().print("\n{s} OPTIMIZED SYSTEM:\n", .{name}); try system.print(std.io.getStdOut().writer());}
    }

    const hf = try hfm.runFull(T, opt.hartree_fock, system, print, allocator); var mp = try mpPost(T, opt, system.*, hf, print, allocator);

    if (opt.gradient != null) mp.G = try edf.gradient(T, opt, system.*, mpFull, name, true, allocator);

    if (print) {
        if (mp.G != null) {try std.io.getStdOut().writer().print("\n{s} GRADIENT:\n", .{name}); try mp.G.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null) mp.H = try edf.hessian(T, opt, system.*, mpFull, name, true, allocator);

    if (print and opt.hessian != null and opt.print.hessian) {
        if (mp.H != null) {try std.io.getStdOut().writer().print("\n{s} HESSIAN:\n", .{name}); try mp.H.?.print(std.io.getStdOut().writer());}
    }

    if (opt.hessian != null and opt.hessian.?.freq) mp.freqs = try prp.freq(T, system.*, mp.H.?, allocator);

    if (print) {
        if (mp.freqs != null) {try std.io.getStdOut().writer().print("\n{s} HARMONIC FREQUENCIES:\n", .{name}); try mp.freqs.?.matrix().print(std.io.getStdOut().writer());}
    }

    if (opt.gradient != null and opt.write.gradient != null) try mp.G.?.write(opt.write.gradient.?);
    if (opt.hessian  != null and opt.write.hessian  != null) try mp.H.?.write(opt.write.hessian.? );

    return mp;
}

/// Function to actually run the Moller-Plesset energy calculation on the provided system without additional calculation.
pub fn mpFull(comptime T: type, opt: inp.MollerPlessetOptions(T), system: System(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    const hf = try hfm.hfFull(T, opt.hartree_fock, system, print, allocator); return mpPost(T, opt, system, hf, print, allocator);
}

/// Function to run the Moller-Plesset energy calculation on the provided system with Hartree-Fock output.
pub fn mpPost(comptime T: type, opt: inp.MollerPlessetOptions(T), system: System(T), hf: out.HartreeFockOutput(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    const nbf = if (opt.hartree_fock.generalized) hf.S_A.rows else 2 * hf.S_A.rows; const nocc = system.getElectrons();

    var F_MS   = try Matrix(T).init(nbf, nbf,                      allocator); defer   F_MS.deinit();
    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();

    try transform(T, &F_MS, &J_MS_A, hf.F_A, hf.J_A, hf.C_A, allocator);

    var E: T = 0;

    if (opt.order >= 2) E += mp2(T, F_MS, J_MS_A, nocc);

    if (print) try std.io.getStdOut().writer().print("\nMP{d} ENERGY: {d:.14}\n", .{opt.order, hf.E + E});

    return .{
        .hf = hf, .E = hf.E + E, .G = null, .H = null, .freqs = null,
    };
}

/// Check for errors in the Moller-Plesset options.
pub fn checkErrors(comptime T: type, opt: inp.MollerPlessetOptions(T)) !void {
    const readInt = opt.hartree_fock.integral.overlap != null or opt.hartree_fock.integral.kinetic != null or opt.hartree_fock.integral.nuclear != null or opt.hartree_fock.integral.coulomb != null;

    if (opt.gradient != null and readInt) return error.CantUseGradientWithProvidedIntegrals;
    if (opt.hessian  != null and readInt) return  error.CantUseHessianWithProvidedIntegrals;

    if (opt.hartree_fock.optimize != null) return error.CantOptimizeOnHartreeFockAndUseMollerPlesset;

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
pub fn transform(comptime T: type, F_MS: *Matrix(T), J_MS_A: *Tensor(T), F_A: Matrix(T), J_A: Tensor(T), C_A: Matrix(T), allocator: std.mem.Allocator) !void {
    if (F_A.rows != F_MS.rows) {try tns.oneAO2MS(T, F_MS,   F_A, C_A           );} else {try tns.oneAO2MO(T, F_MS,   F_A, C_A           );}
    if (F_A.rows != F_MS.rows) {try tns.twoAO2MS(T, J_MS_A, J_A, C_A, allocator);} else {try tns.twoAO2MO(T, J_MS_A, J_A, C_A, allocator);}

    for (0..J_MS_A.shape[0]) |i| for (0..J_MS_A.shape[1]) |j| for (0..J_MS_A.shape[2]) |k| for (0..J_MS_A.shape[3]) |l| {
        J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS_A.at(&[_]usize{i, j, k, l}) - J_MS_A.at(&[_]usize{i, l, k, j});
    };
}

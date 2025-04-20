//! Module to perform the Moller-Plesset calculations.

const std = @import("std");

const inp = @import("input.zig"      );
const hfm = @import("hartreefock.zig");
const out = @import("output.zig"     );
const ten = @import("tensor.zig"     );
const tns = @import("transform.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// Main function to run the Moller-Plesset calculations.
pub fn run(comptime T: type, opt: inp.MollerPlessetOptions(T), print: bool, allocator: std.mem.Allocator) !out.MollerPlessetOutput(T) {
    if (opt.order != 2) return error.PerturbationOrderNotImplemented;

    const hf = try hfm.run(T, opt.hartree_fock, print, allocator); const nbf = 2 * hf.S_AO.rows; const nocc = 2 * hf.system.nocc;

    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();
    
    try transform(T, &J_MS_A, hf.J_AO, hf.C_MO, allocator);

    var E: T = 0;

    if (opt.order >= 2) E += mp2(T, hf.E_MO, J_MS_A, nocc);

    if (print) try std.io.getStdOut().writer().print("\nMP{d} ENERGY: {d:.14}\n", .{opt.order, hf.E + E});

    return .{
        .E = hf.E + E,
    };
}

/// Returns the second-order Moller-Plesset energy.
pub fn mp2(comptime T: type, E_MO: Matrix(T), J_MS_A: Tensor(T), nocc: usize) T {
    var E: T = 0;

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..J_MS_A.shape[0]) |a| for (nocc..J_MS_A.shape[0]) |b| {
        E += 0.25 * J_MS_A.at(&[_]usize{i, j, a, b}) * J_MS_A.at(&[_]usize{a, b, i, j}) / (E_MO.at(i / 2, i / 2) + E_MO.at(j / 2, j / 2) - E_MO.at(a / 2, a / 2) - E_MO.at(b / 2, b / 2));
    };

    return E;
}

/// Function to perform all integrals transformations used in the Moller-Plesset calculations.
pub fn transform(comptime T: type, J_MS_A: *Tensor(T), J_AO: Tensor(T), C_MO: Matrix(T), allocator: std.mem.Allocator) !void {
    var J_AS = try Tensor(T).init(&[_]usize{J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0]}, allocator); defer J_AS.deinit();
    var J_MS = try Tensor(T).init(&[_]usize{J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0]}, allocator); defer J_MS.deinit();

    var C_MS = try Matrix(T).init(J_MS_A.shape[0], J_MS_A.shape[0], allocator); defer C_MS.deinit();

    tns.cfsMO2MS(T, &C_MS, C_MO); tns.twoAO2AS(T, &J_AS, J_AO); try tns.twoAO2MO(T, &J_MS, &J_AS, C_MS);

    for (0..J_MS.shape[0]) |i| for (0..J_MS.shape[1]) |j| for (0..J_MS.shape[2]) |k| for (0..J_MS.shape[3]) |l| {
        J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS.at(&[_]usize{i, j, k, l}) - J_MS.at(&[_]usize{i, l, k, j});
    };
}

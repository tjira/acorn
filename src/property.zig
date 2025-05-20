//! File to calculate properties from a quantum mechanical calculation

const std = @import("std");

const AN2M = @import("constant.zig").AN2M;

const cnt = @import("constant.zig");
const lpk = @import("lapack.zig"  );
const mth = @import("math.zig"    );

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Vector = @import("vector.zig").Vector;

/// Returns the harmonic frequencies of the system.
pub fn freq(comptime T: type, system: System(T), H: Matrix(T), allocator: std.mem.Allocator) !Vector(T) {
    const freqs = try Vector(T).init(H.rows, allocator); freqs.fill(0); var HM = try H.clone(); defer HM.deinit();

    for (system.atoms.data) |atom| {if (@as(usize, @intFromFloat(atom)) - 1 > AN2M.len) return error.UnknownMassForAtom;}

    for (0..H.rows) |i| for (0..H.rows) |j| {

        const an1 = @as(usize, @intFromFloat(system.atoms.at(i / 3)));
        const an2 = @as(usize, @intFromFloat(system.atoms.at(j / 3)));

        HM.ptr(i, j).* /= std.math.sqrt(AN2M[an1 - 1] * AN2M[an2 - 1]);
    };

    var HMJ = try Matrix(T).init(HM.rows, HM.cols, allocator); defer HMJ.deinit();
    var HMC = try Matrix(T).init(HM.rows, HM.cols, allocator); defer HMC.deinit();

    lpk.dsyevd(&HMJ, &HMC, HM);

    const factor = 5e-3 / std.math.pi / cnt.c * std.math.sqrt(cnt.Eh / cnt.amu / cnt.a0 / cnt.a0);

    for (0..freqs.rows) |i| freqs.ptr(i).* = factor * mth.sgn(HMJ.at(i, i)) * std.math.sqrt(@abs(HMJ.at(i, i)));

    return freqs;
}

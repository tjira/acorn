//! Functions to generate CAS determinants.

const std = @import("std");

const mth = @import("math.zig");

const Matrix = @import("matrix.zig").Matrix;

/// Generates the determinants for the CI calculations including singlet and triplet excitations.
pub fn generateAllDeterminants(nbf: usize, nel: usize, allocator: std.mem.Allocator) !Matrix(usize) {
    const D = try Matrix(usize).init(mth.comb(nbf, nel), nel, allocator);

    for (0..D.cols) |i| D.ptr(0, i).* = i;

    for (1..D.rows) |i| {

        @memcpy(D.row(i).data, D.row(i - 1).data); var index: usize = undefined;

        for (0..D.cols) |j| if (D.at(i, D.cols - j - 1) != nbf - j - 1) {
            D.ptr(i, D.cols - j - 1).* += 1; index = D.cols - j - 1; break;
        };

        for (index + 1..D.cols) |j| D.ptr(i, j).* = D.at(i, j - 1) + 1;
    }

    return D;
}

/// Generates CAS determinants for the CI calculations.
pub fn generateCasDeterminants(nel: usize, nbf: usize, nocc: usize, spin: bool, allocator: std.mem.Allocator) !Matrix(usize) {
    const CASD = if (spin) try generateAllDeterminants(nbf, nel, allocator) else try generateAllSingletDeterminants(nbf, nel, allocator); defer CASD.deinit();

    var D = try Matrix(usize).init(CASD.rows, nocc, allocator);

    for (0..CASD.rows) |i| {
        for (0..nocc - nel) |j| {D.ptr(i, j).* = j;} for (0..CASD.cols) |j| D.ptr(i, j + nocc - nel).* = CASD.at(i, j) + nocc - nel;
    }

    return D;
}

/// Generates the determinants for the CI calculations including only singlet excitations.
pub fn generateAllSingletDeterminants(nbf: usize, nel: usize, allocator: std.mem.Allocator) !Matrix(usize) {
    const data_alpha = try allocator.alloc(usize, nbf / 2); defer allocator.free(data_alpha);
    const data_beta  = try allocator.alloc(usize, nbf / 2); defer allocator.free(data_beta );
    
    for (0..data_alpha.len) |i| data_alpha[i] = 2 * i + 0;
    for (0..data_beta.len ) |i| data_beta[i]  = 2 * i + 1;

    const dets_alpha = try mth.combinations(usize, data_alpha, nel / 2, allocator); defer dets_alpha.deinit();
    const dets_beta  = try mth.combinations(usize, data_beta,  nel / 2, allocator); defer  dets_beta.deinit();

    const D = try Matrix(usize).init(dets_alpha.items.len * dets_beta.items.len, nel, allocator);

    for (0..dets_alpha.items.len) |i| for (0..dets_beta.items.len) |j| for (0..nel / 2) |k| {
        D.ptr(i * dets_beta.items.len + j, k           ).* = dets_alpha.items[i][k];
        D.ptr(i * dets_beta.items.len + j, k + nel / 2).* =  dets_beta.items[j][k];
    };

    return D;
}

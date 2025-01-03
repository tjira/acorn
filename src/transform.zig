//! Module for the integral transformation.

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// Transforms coefficients obtained from the HF calculation from the MO basis to the MS basis.
pub fn cfsMO2MS(comptime T: type, C_MS: *Matrix(T), C_MO: Matrix(T)) void {
    C_MS.fill(0);

    for (0..C_MO.rows) |i| for (0..C_MO.cols) |j| {
        C_MS.ptr(i            , 2 * j    ).* = C_MO.at(i, j);
        C_MS.ptr(i + C_MO.rows, 2 * j + 1).* = C_MO.at(i, j);
    };
}

/// Transforms the one-electron integrals from the AO basis to the AS basis.
pub fn oneAO2AS(comptime T: type, A_AS: *Matrix(T), A_AO: Matrix(T)) void {
    A_AS.fill(0);

    for (0..A_AO.rows) |i| for (0..A_AO.cols) |j| {
        A_AS.ptr(i            , j            ).* = A_AO.at(i, j);
        A_AS.ptr(i + A_AO.rows, j + A_AO.cols).* = A_AO.at(i, j);
    };
}

/// Transforms the one-electron integrals from the AO basis to the MO basis or from the AS basis to the MS basis.
pub fn oneAO2MO(comptime T: type, A_MO: *Matrix(T), A_AO: Matrix(T), C_MO: Matrix(T)) void {
    A_MO.fill(0);

    for (0..A_AO.rows) |i| for (0..A_AO.cols) |j| for (0..A_MO.rows) |p| for (0..A_MO.cols) |q| {
        A_MO.ptr(p, q).* += C_MO.at(i, p) * A_AO.at(i, j) * C_MO.at(j, q);
    };
}

/// Transforms the two-electron integrals from the AO basis to the AS basis.
pub fn twoAO2AS(comptime T: type, A_AS: *Tensor(T), A_AO: Tensor(T)) void {
    A_AS.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| {
        A_AS.ptr(&[_]usize{i,                 j,                 k                , l                }).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i,                 j                , k + A_AO.shape[2], l + A_AO.shape[3]}).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i + A_AO.shape[0], j + A_AO.shape[1], k,                 l                }).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i + A_AO.shape[0], j + A_AO.shape[1], k + A_AO.shape[2], l + A_AO.shape[3]}).* = A_AO.at(&[_]usize{i, j, k, l});
    };
}

/// Transforms the two-electron integrals from the AO basis to the MO basis or from the AS basis to the MS basis.
pub fn twoAO2MO(comptime T: type, A_MO: *Tensor(T), A_AO: *Tensor(T), C_MO: Matrix(T)) void {
    @memcpy(A_MO.data, A_AO.data); A_AO.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| for (0..A_MO.shape[0]) |p| {
        A_AO.ptr(&[_]usize{p, j, k, l}).* += A_MO.at(&[_]usize{i, j, k, l}) * C_MO.at(i, p);
    };

    A_MO.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| for (0..A_MO.shape[1]) |q| {
        A_MO.ptr(&[_]usize{i, q, k, l}).* += A_AO.at(&[_]usize{i, j, k, l}) * C_MO.at(j, q);
    };

    A_AO.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| for (0..A_MO.shape[2]) |r| {
        A_AO.ptr(&[_]usize{i, j, r, l}).* += A_MO.at(&[_]usize{i, j, k, l}) * C_MO.at(k, r);
    };

    A_MO.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| for (0..A_MO.shape[3]) |s| {
        A_MO.ptr(&[_]usize{i, j, k, s}).* += A_AO.at(&[_]usize{i, j, k, l}) * C_MO.at(l, s);
    };
}


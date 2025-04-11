//! This file contains wrappers for the LAPACK library.

const lapack = @cImport(@cInclude("lapacke.h"));

const std = @import("std");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// Solve the linear system Ax = b using the LU decomposition of A. The result is stored in x. The LU decomposition is stored in ALU and the pivot indices are stored in p.
pub fn dgesv(x: *Vector(f64), ALU: *Matrix(f64), p: *Vector(i32), A: Matrix(f64), b: Vector(f64)) void {
    const n: i32 = @intCast(A.rows); @memcpy(ALU.data, A.data); @memcpy(x.data, b.data);

    _ = lapack.LAPACKE_dgesv(lapack.LAPACK_COL_MAJOR, n, 1, &ALU.data[0], n, &p.data[0], &x.data[0], n);
}

/// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
pub fn dsyevd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64)) void {
    const n: i32 = @intCast(A.rows); @memcpy(C.data, A.data); J.fill(0);

    _ = lapack.LAPACKE_dsyevd(lapack.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));
}

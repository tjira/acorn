//! This file contains wrappers for the LAPACK library.

const lapacke = @cImport(@cInclude("lapacke.h"));

const std = @import("std");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// Calculate the condition number of a matrix A. The result returned.
pub fn dgecon(ALU: Matrix(f64), onorm: f64) f64 {
    const n: i32 = @intCast(ALU.rows); var rcond: f64 = 0;

    _ = lapacke.LAPACKE_dgecon(lapacke.LAPACK_ROW_MAJOR, '1', n, &ALU.data[0], n, onorm, &rcond);

    return rcond;
}

/// Calculate the Schur decomposition of matrix A. The JR and JI will contain real and imaginary parts of calculated eigenvalues.
pub fn dgees(Q: *Matrix(f64), D: *Matrix(f64), A: Matrix(f64), JR: *Vector(f64), JI: *Vector(f64)) void {
    const n: i32 = @intCast(A.rows); var sdim: i32 = undefined; @memcpy(D.data, A.data);

    _ = lapacke.LAPACKE_dgees(lapacke.LAPACK_ROW_MAJOR, 'V', 'N', null, n, &D.data[0], n, &sdim, &JR.data[0], &JI.data[0], &Q.data[0], n);
}

/// Solve the linear system Ax = b using the LU decomposition of A. The result is stored in x. The LU decomposition is stored in ALU and the pivot indices are stored in p.
pub fn dgesv(x: *Vector(f64), ALU: *Matrix(f64), p: *Vector(i32), A: Matrix(f64), b: Vector(f64)) void {
    const n: i32 = @intCast(A.rows); @memcpy(ALU.data, A.data); @memcpy(x.data, b.data);

    _ = lapacke.LAPACKE_dgesv(lapacke.LAPACK_ROW_MAJOR, n, 1, &ALU.data[0], n, &p.data[0], &x.data[0], 1);
}

/// Calculate the LU decomposition of a matrix A. The result is stored in ALU and the pivot indices are stored in p.
pub fn dgetrf(ALU: *Matrix(f64), p: *Vector(i32), A: Matrix(f64)) void {
    const m: i32 = @intCast(A.rows); const n: i32 = @intCast(A.cols); @memcpy(ALU.data, A.data);

    _ = lapacke.LAPACKE_dgetrf(lapacke.LAPACK_ROW_MAJOR, m, n, &ALU.data[0], m, &p.data[0]);
}

/// Calculate the inverse from the LU decomposition of a matrix A. The result is stored in the AIN matrix.
pub fn dgetri(AIN: *Matrix(f64), ALU: Matrix(f64), p: Vector(i32)) void {
    const n: i32 = @intCast(ALU.rows); @memcpy(AIN.data, ALU.data);

    _ = lapacke.LAPACKE_dgetri(lapacke.LAPACK_ROW_MAJOR, n, &AIN.data[0], n, &p.data[0]);
}

/// Calculates the specified norm of a matrix A. The result is returned.
pub fn dlange(A: Matrix(f64), mode: u8) f64 {
    const m: i32 = @intCast(A.rows); const n: i32 = @intCast(A.cols);

    return lapacke.LAPACKE_dlange(lapacke.LAPACK_ROW_MAJOR, mode, m, n, &A.data[0], m);
}

/// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
pub fn dsyevd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64)) void {
    const n: i32 = @intCast(A.rows); @memcpy(C.data, A.data); J.fill(0);

    _ = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));
}

/// Finds the eigenvalues and eigenvectors of a symmetric-definite generalized eigenproblem A*x = λ*B*x. The eigenvalues are stored in J and the eigenvectors are stored in C. The upper triangular part of Cholesky decmposition will be stored in CH.
pub fn dsygvd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64), B: Matrix(f64), CH: *Matrix(f64)) void {
    const n: i32 = @intCast(A.rows); @memcpy(C.data, A.data); @memcpy(CH.data, B.data); J.fill(0);

    _ = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, &C.data[0], n, &CH.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));
}

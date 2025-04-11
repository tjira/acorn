//! This file contains wrappers for the BLAS library.

const std = @import("std"); const Complex = std.math.Complex;

const blas = @cImport(@cInclude("cblas.h"));

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AH and BH flags control whether the matrices A and B are Hermitian transposed or not.
pub fn zgemm(C: *Matrix(Complex(f64)), A: Matrix(Complex(f64)), AH: bool, B: Matrix(Complex(f64)), BH: bool) void {
    const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

    const ah: c_uint = if (AH) blas.CblasConjTrans else blas.CblasNoTrans;
    const bh: c_uint = if (BH) blas.CblasConjTrans else blas.CblasNoTrans;

    _ = blas.cblas_zgemm(blas.CblasRowMajor, ah, bh, m, n, k, &Complex(f64).init(1.0, 0.0), &A.data[0], m, &B.data[0], n, &Complex(f64).init(0.0, 0.0), &C.data[0], n);
}

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AT and BT flags control whether the matrices A and B are transposed or not.
pub fn dgemm(C: *Matrix(f64), A: Matrix(f64), AT: bool, B: Matrix(f64), BT: bool) void {
    const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

    const at: c_uint = if (AT) blas.CblasTrans else blas.CblasNoTrans;
    const bt: c_uint = if (BT) blas.CblasTrans else blas.CblasNoTrans;

    _ = blas.cblas_dgemm(blas.CblasRowMajor, at, bt, m, n, k, 1.0, &A.data[0], m, &B.data[0], n, 0.0, &C.data[0], n);
}

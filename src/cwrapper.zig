//! Wrappers for the C/C++ libraries.

const cblas         = @cImport(@cInclude("cblas.h"            ));
const eigen         = @cImport(@cInclude("eigen.h"            ));
const exprtk        = @cImport(@cInclude("exprtk.h"           ));
const fftw          = @cImport(@cInclude("fftw3.h"            ));
const gsl_sf_fact   = @cImport(@cInclude("gsl/gsl_sf_fact.h"  ));
const gsl_sf_gamma  = @cImport(@cInclude("gsl/gsl_sf_gamma.h" ));
const gsl_sf_hyperg = @cImport(@cInclude("gsl/gsl_sf_hyperg.h"));
const lapacke       = @cImport(@cInclude("lapacke.h"          ));
const libint        = @cImport(@cInclude("libint.h"           ));

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Compile the algebrainc expression string with the given number of variables. Variables must be in the form of r1, r2, ..., rn.
pub fn compileExpression(string: []const u8, nvars: usize, allocator: std.mem.Allocator) !*anyopaque {
    const buffer = try allocator.dupeZ(u8, string); defer allocator.free(buffer); return exprtk.compile(buffer, nvars).?;
}

pub fn contract(C: anytype, A: anytype, B: anytype, pairs: []const i32) !void {
    var dC: []const usize = undefined; var rC: usize = undefined;
    var dA: []const usize = undefined; var rA: usize = undefined;
    var dB: []const usize = undefined; var rB: usize = undefined;

    if (comptime std.meta.fieldIndex(@TypeOf(C.*), "shape") != null) {dC = C.shape; rC = C.shape.len;} else {dC = &[_]usize{C.rows, C.cols}; rC = 2;}
    if (comptime std.meta.fieldIndex(@TypeOf(A  ), "shape") != null) {dA = A.shape; rA = A.shape.len;} else {dA = &[_]usize{A.rows, A.cols}; rA = 2;}
    if (comptime std.meta.fieldIndex(@TypeOf(B  ), "shape") != null) {dB = B.shape; rB = B.shape.len;} else {dB = &[_]usize{B.rows, B.cols}; rB = 2;}

    return eigen.contract(&C.data[0], &dC[0], rC, &A.data[0], &dA[0], rA, &B.data[0], &dB[0], rB, &pairs[0], pairs.len / 2);
}

/// Calculate the Coulomb integrals.
pub fn coulombIntegrals(ints: []f64, system: System(f64), basis: std.ArrayList(f64)) void {
    libint.coulomb(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
}

/// Deinitialize the compiled expression.
pub fn deinitExpression(expr: *anyopaque) void {
    exprtk.deinit(expr);
}

/// Double factorial function on unsigned integers.
pub fn dfact(n: u32) f64 {
    return gsl_sf_gamma.gsl_sf_doublefact(n);
}

/// Calculate the condition number of a matrix A. The result returned.
pub fn dgecon(ALU: Matrix(f64), onorm: f64) f64 {
    const n: i32 = @intCast(ALU.rows); var rcond: f64 = 0;

    _ = lapacke.LAPACKE_dgecon(lapacke.LAPACK_ROW_MAJOR, '1', n, &ALU.data[0], n, onorm, &rcond);

    return rcond;
}

/// Calculate the Schur decomposition of matrix A. The JR and JI will contain real and imaginary parts of calculated eigenvalues.
pub fn dgees(Q: *Matrix(f64), D: *Matrix(f64), A: Matrix(f64), JR: *Vector(f64), JI: *Vector(f64)) void {
    const n: i32 = @intCast(A.rows); var sdim: i32 = undefined; A.memcpy(D.*);

    _ = lapacke.LAPACKE_dgees(lapacke.LAPACK_ROW_MAJOR, 'V', 'N', null, n, &D.data[0], n, &sdim, &JR.data[0], &JI.data[0], &Q.data[0], n);
}

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AT and BT flags control whether the matrices A and B are transposed or not.
pub fn dgemm(C: *Matrix(f64), A: Matrix(f64), AT: bool, B: Matrix(f64), BT: bool) void {
    const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

    const at: c_uint = if (AT) cblas.CblasTrans else cblas.CblasNoTrans;
    const bt: c_uint = if (BT) cblas.CblasTrans else cblas.CblasNoTrans;

    _ = cblas.cblas_dgemm(cblas.CblasRowMajor, at, bt, m, n, k, 1.0, &A.data[0], m, &B.data[0], n, 0.0, &C.data[0], n);
}

/// Solve the linear system Ax = b using the LU decomposition of A. The result is stored in x. The LU decomposition is stored in ALU and the pivot indices are stored in p.
pub fn dgesv(x: *Vector(f64), ALU: *Matrix(f64), p: *Vector(i32), A: Matrix(f64), b: Vector(f64)) void {
    const n: i32 = @intCast(A.rows); A.memcpy(ALU.*); b.memcpy(x.*);

    _ = lapacke.LAPACKE_dgesv(lapacke.LAPACK_ROW_MAJOR, n, 1, &ALU.data[0], n, &p.data[0], &x.data[0], 1);
}

/// Calculate the LU decomposition of a matrix A. The result is stored in ALU and the pivot indices are stored in p.
pub fn dgetrf(ALU: *Matrix(f64), p: *Vector(i32), A: Matrix(f64)) void {
    const m: i32 = @intCast(A.rows); const n: i32 = @intCast(A.cols); A.memcpy(ALU.*);

    _ = lapacke.LAPACKE_dgetrf(lapacke.LAPACK_ROW_MAJOR, m, n, &ALU.data[0], m, &p.data[0]);
}

/// Calculate the inverse from the LU decomposition of a matrix A. The result is stored in the AIN matrix.
pub fn dgetri(AIN: *Matrix(f64), ALU: Matrix(f64), p: Vector(i32)) void {
    const n: i32 = @intCast(ALU.rows); ALU.memcpy(AIN.*);

    _ = lapacke.LAPACKE_dgetri(lapacke.LAPACK_ROW_MAJOR, n, &AIN.data[0], n, &p.data[0]);
}

/// Calculates the specified norm of a matrix A. The result is returned.
pub fn dlange(A: Matrix(f64), mode: u8) f64 {
    const m: i32 = @intCast(A.rows); const n: i32 = @intCast(A.cols);

    return lapacke.LAPACKE_dlange(lapacke.LAPACK_ROW_MAJOR, mode, m, n, &A.data[0], m);
}

/// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
pub fn dsyevd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64)) void {
    const n: i32 = @intCast(A.rows); A.memcpy(C.*); J.fill(0);

    _ = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));
}

/// Finds the eigenvalues and eigenvectors of a symmetric-definite generalized eigenproblem A*x = λ*B*x. The eigenvalues are stored in J and the eigenvectors are stored in C. The upper triangular part of Cholesky decmposition will be stored in CH.
pub fn dsygvd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64), B: Matrix(f64), CH: *Matrix(f64)) void {
    const n: i32 = @intCast(A.rows); A.memcpy(C.*); B.memcpy(CH.*); J.fill(0);

    _ = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, &C.data[0], n, &CH.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));
}

/// Evaluate the compiled expression at the specified point.
pub fn evaluateExpression(expr: *anyopaque, vars: Vector(f64)) f64 {
    return exprtk.evaluate(expr, vars.data.ptr);
}

pub fn fftwnd(in: []std.math.Complex(f64), shape: []const u32, factor: i32) void {
    const plan = fftw.fftw_plan_dft(@intCast(shape.len), &@as([]const i32, @ptrCast(shape))[0], &@as([][2]f64, @ptrCast(in))[0], &@as([][2]f64, @ptrCast(in))[0], factor, fftw.FFTW_ESTIMATE);

    fftw.fftw_execute(plan); fftw.fftw_destroy_plan(plan);

    if (factor > 0) for (0..in.len) |i| {
        in[i] = in[i].div(std.math.Complex(f64).init(@as(f64, @floatFromInt(in.len)), 0));
    };
}

/// Gamma function.
pub fn gamma(x: f64) f64 {
    return gsl_sf_gamma.gsl_sf_gamma(x);
}

/// AArgument of Complex Gamma function.
pub fn gammaArg(z: std.math.Complex(f64)) f64 {
    var ln  = gsl_sf_gamma.gsl_sf_result{};
    var arg = gsl_sf_gamma.gsl_sf_result{};

    _ = gsl_sf_gamma.gsl_sf_lngamma_complex_e(z.re, z.im, &ln, &arg);

    return arg.val;
}

/// Regularized lower incomplete gamma function.
pub fn gammainc(a: f64, x: f64) f64 {
    return gsl_sf_gamma.gsl_sf_gamma_inc_P(a, x);
}

/// Kummer confluent hypergeometric function.
pub fn hyp1f1(a: f64, b: f64, x: f64) f64 {
    return gsl_sf_hyperg.gsl_sf_hyperg_1F1(a, b, x);
}

/// Calculate the kinetic integrals.
pub fn kineticIntegrals(ints: []f64, system: System(f64), basis: std.ArrayList(f64)) void {
    libint.kinetic(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
}

/// Logarithm of a matrix.
pub fn logm(B: *Matrix(f64), A: Matrix(f64)) void {
    eigen.logm(&B.data[0], &A.data[0], A.rows);
}

/// Calculate the nuclear integrals.
pub fn nuclearIntegrals(ints: []f64, system: System(f64), basis: std.ArrayList(f64)) void {
    libint.nuclear(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
}

/// Calculate the overlap integrals.
pub fn overlapIntegrals(ints: []f64, system: System(f64), basis: std.ArrayList(f64)) void {
    libint.overlap(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
}

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AH and BH flags control whether the matrices A and B are Hermitian transposed or not.
pub fn zgemm(C: *Matrix(std.math.Complex(f64)), A: Matrix(std.math.Complex(f64)), AH: bool, B: Matrix(std.math.Complex(f64)), BH: bool) void {
    const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

    const ah: c_uint = if (AH) cblas.CblasConjTrans else cblas.CblasNoTrans;
    const bh: c_uint = if (BH) cblas.CblasConjTrans else cblas.CblasNoTrans;

    _ = cblas.cblas_zgemm(cblas.CblasRowMajor, ah, bh, m, n, k, &std.math.Complex(f64).init(1.0, 0.0), &A.data[0], m, &B.data[0], n, &std.math.Complex(f64).init(0.0, 0.0), &C.data[0], n);
}

/// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
pub fn zheevd(J: *Matrix(std.math.Complex(f64)), C: *Matrix(std.math.Complex(f64)), A: Matrix(std.math.Complex(f64))) !void {
    A.memcpy(C.*); J.fill(std.math.Complex(f64).init(0, 0));

    const AR = try Matrix(f64).init(2 * A.rows, 2 * A.cols, A.allocator);
    var   JR = try Matrix(f64).init(2 * J.rows, 2 * J.cols, J.allocator);
    var   CR = try Matrix(f64).init(2 * C.rows, 2 * C.cols, C.allocator);

    for (0..A.rows) |i| for (0..A.cols) |j| {
        AR.ptr(i,          j         ).* =  A.at(i, j).re;
        AR.ptr(i,          j + A.cols).* = -A.at(i, j).im;
        AR.ptr(i + A.rows, j         ).* =  A.at(i, j).im;
        AR.ptr(i + A.rows, j + A.cols).* =  A.at(i, j).re;
    };

    dsyevd(&JR, &CR, AR);

    for (0..J.rows) |i| J.ptr(i, i).*.re = JR.at(2 * i, 2 * i);

    for (0..C.rows) |i| for (0..C.cols) |j| {
        C.ptr(i, j).*.re = CR.at(i,          2 * j);
        C.ptr(i, j).*.im = CR.at(i + C.rows, 2 * j);
    };
}

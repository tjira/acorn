//! Wrappers for the C/C++ libraries.

const cblas        = @cImport(@cInclude("cblas.h"           ));
const eigen        = @cImport(@cInclude("eigen.h"           ));
const exprtk       = @cImport(@cInclude("exprtk.h"          ));
const fftw         = @cImport(@cInclude("fftw3.h"           ));
const gsl_sf_gamma = @cImport(@cInclude("gsl/gsl_sf_gamma.h"));
const lapacke      = @cImport(@cInclude("lapacke.h"         ));
const libint       = @cImport(@cInclude("libint.h"          ));

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Expression class.
pub fn Expression(comptime T: type) type {
    return struct {
        expression: *anyopaque,

        /// Initializer function.
        pub fn init(string: []const u8, nvars: usize, allocator: std.mem.Allocator) !Expression(T) {
            const buffer = try allocator.dupeZ(u8, string); defer allocator.free(buffer); return .{.expression = exprtk.compile(buffer, nvars).?};
        }

        /// Deinitialize the expression.
        pub fn deinit(self: Expression(T)) void {
            exprtk.deinit(self.expression);
        }

        /// Evaluate the compiled expression at the specified point.
        pub fn evaluate(self: Expression(T), vars: Vector(T), time: T) T {
            return exprtk.evaluate(self.expression, vars.data.ptr, time);
        }
    };
}

/// Libint struct.
pub fn Libint(comptime T: type) type {
    return struct {

        /// Calculate the Coulomb integrals.
        pub fn coulomb(ints: []T, system: System(T), basis: std.ArrayList(T)) void {
            libint.coulomb(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the kinetic integrals.
        pub fn kinetic(ints: []T, system: System(T), basis: std.ArrayList(T)) void {
            libint.kinetic(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the nuclear integrals.
        pub fn nuclear(ints: []T, system: System(T), basis: std.ArrayList(T)) void {
            libint.nuclear(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the overlap integrals.
        pub fn overlap(ints: []T, system: System(T), basis: std.ArrayList(T)) void {
            libint.overlap(&ints[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }
    };
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

/// Calculate the condition number of a matrix A, provided its LU decomposition and 1-norm. The result returned.
pub fn dgecon(ALU: Matrix(f64), onorm: f64) !f64 {
    const n: i32 = @intCast(ALU.rows); var rcond: f64 = 0;

    const info = lapacke.LAPACKE_dgecon(lapacke.LAPACK_ROW_MAJOR, '1', n, &ALU.data[0], n, onorm, &rcond);

    return if (info == 0) rcond else error.ErrorInConditionNumberCalculation;
}

/// Calculate the Schur decomposition of matrix A. The JR and JI will contain real and imaginary parts of calculated eigenvalues.
pub fn dgees(T: *Matrix(f64), Q: *Matrix(f64), A: Matrix(f64), JR: *Matrix(f64), JI: *Matrix(f64)) !void {
    const n: i32 = @intCast(A.rows); var sdim: i32 = undefined; A.memcpy(T.*);

    const info = lapacke.LAPACKE_dgees(lapacke.LAPACK_ROW_MAJOR, 'V', 'N', null, n, &T.data[0], n, &sdim, &JR.data[0], &JI.data[0], &Q.data[0], n);

    for (0..JR.rows) |i| std.mem.swap(f64, JR.ptr(i, i), JR.ptr(0, i));
    for (0..JI.rows) |i| std.mem.swap(f64, JI.ptr(i, i), JI.ptr(0, i));

    if (info != 0) return error.ErrorInRealSchurDecomposition;
}

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AT and BT flags control whether the matrices A and B are transposed or not.
pub fn dgemm(C: *Matrix(f64), A: Matrix(f64), AT: bool, B: Matrix(f64), BT: bool) void {
    const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

    const at: c_uint = if (AT) cblas.CblasTrans else cblas.CblasNoTrans;
    const bt: c_uint = if (BT) cblas.CblasTrans else cblas.CblasNoTrans;

    cblas.cblas_dgemm(cblas.CblasRowMajor, at, bt, m, n, k, 1.0, &A.data[0], m, &B.data[0], n, 0.0, &C.data[0], n);
}

/// Solve the linear system Ax = b using the LU decomposition of A. The result is stored in x. The LU decomposition is stored in ALU and the pivot indices are stored in p.
pub fn dgesv(x: *Vector(f64), ALU: *Matrix(f64), p: *Vector(i32), A: Matrix(f64), b: Vector(f64)) !void {
    const n: i32 = @intCast(A.rows); A.memcpy(ALU.*); b.memcpy(x.*);

    const info = lapacke.LAPACKE_dgesv(lapacke.LAPACK_ROW_MAJOR, n, 1, &ALU.data[0], n, &p.data[0], &x.data[0], 1);

    if (info != 0) return error.ErrorInLinearSystemSolution;
}

/// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
pub fn dsyevd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64)) !void {
    const n: i32 = @intCast(A.rows); A.memcpy(C.*); J.fill(0);

    const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));

    if (info != 0) return error.ErrorInEigenvalueCalculation;
}

/// Finds the eigenvalues and eigenvectors of a symmetric-definite generalized eigenproblem A*x = λ*B*x. The eigenvalues are stored in J and the eigenvectors are stored in C. The upper triangular part of Cholesky decmposition will be stored in CH.
pub fn dsygvd(J: *Matrix(f64), C: *Matrix(f64), A: Matrix(f64), B: Matrix(f64), CH: *Matrix(f64)) !void {
    const n: i32 = @intCast(A.rows); A.memcpy(C.*); B.memcpy(CH.*); J.fill(0);

    const info = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, &C.data[0], n, &CH.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(f64, J.ptr(i, i), J.ptr(0, i));

    if (info != 0) return error.ErrorInEigenvalueCalculation;
}

pub fn fftwnd(in: []std.math.Complex(f64), shape: []const u32, factor: i32) void {
    const data = @as([][2]f64, @ptrCast(in));

    const plan = fftw.fftw_plan_dft(@intCast(shape.len), &@as([]const i32, @ptrCast(shape))[0], &data[0], &data[0], factor, fftw.FFTW_ESTIMATE);

    fftw.fftw_execute(plan); fftw.fftw_destroy_plan(plan);

    if (factor > 0) for (0..in.len) |i| {
        in[i] = in[i].div(std.math.Complex(f64).init(@as(f64, @floatFromInt(in.len)), 0));
    };
}

/// Argument of Complex Gamma function.
pub fn gammaArg(z: std.math.Complex(f64)) !f64 {
    var lnr = gsl_sf_gamma.gsl_sf_result{};
    var arg = gsl_sf_gamma.gsl_sf_result{};

    const status = gsl_sf_gamma.gsl_sf_lngamma_complex_e(z.re, z.im, &lnr, &arg);

    return if (status == 0) arg.val else return error.ErrorInGammaFunctionCalculation;
}

/// Logarithm of a matrix.
pub fn logm(B: *Matrix(f64), A: Matrix(f64)) void {
    eigen.logm(&B.data[0], &A.data[0], A.rows);
}

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AH and BH flags control whether the matrices A and B are Hermitian transposed or not.
pub fn zgemm(C: *Matrix(std.math.Complex(f64)), A: Matrix(std.math.Complex(f64)), AH: bool, B: Matrix(std.math.Complex(f64)), BH: bool) void {
    const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

    const ah: c_uint = if (AH) cblas.CblasConjTrans else cblas.CblasNoTrans;
    const bh: c_uint = if (BH) cblas.CblasConjTrans else cblas.CblasNoTrans;

    cblas.cblas_zgemm(cblas.CblasRowMajor, ah, bh, m, n, k, &std.math.Complex(f64).init(1.0, 0.0), &A.data[0], m, &B.data[0], n, &std.math.Complex(f64).init(0.0, 0.0), &C.data[0], n);
}

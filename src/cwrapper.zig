//! Wrappers for the C/C++ libraries.

const cblas        = @cImport(@cInclude("cblas.h"  ));
const eigen        = @cImport(@cInclude("eigen.h"  ));
const exprtk       = @cImport(@cInclude("exprtk.h" ));
const lapacke      = @cImport(@cInclude("lapacke.h"));
const libint       = @cImport(@cInclude("libint.h" ));

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Blas struct.
pub fn Blas(comptime T: type) type {
    return struct {

        /// Calculate the matrix-matrix product C = A * B and store the result in C. The AT and BT flags control whether the matrices A and B are transposed or not.
        pub fn dgemm(C: *Matrix(T), A: Matrix(T), AT: bool, B: Matrix(T), BT: bool) void {
            const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

            const at: c_uint = if (AT) cblas.CblasTrans else cblas.CblasNoTrans;
            const bt: c_uint = if (BT) cblas.CblasTrans else cblas.CblasNoTrans;

            const LDA = if (!AT) k else m;
            const LDB = if (!BT) n else k;

            cblas.cblas_dgemm(cblas.CblasRowMajor, at, bt, m, n, k, 1.0, &A.data[0], LDA, &B.data[0], LDB, 0.0, &C.data[0], n);
        }

        /// Calculate the matrix-matrix product C = A * B and store the result in C. The AH and BH flags control whether the matrices A and B are Hermitian transposed or not.
        pub fn zgemm(C: *Matrix(std.math.Complex(T)), A: Matrix(std.math.Complex(T)), AH: bool, B: Matrix(std.math.Complex(T)), BH: bool) void {
            const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

            const ah: c_uint = if (AH) cblas.CblasConjTrans else cblas.CblasNoTrans;
            const bh: c_uint = if (BH) cblas.CblasConjTrans else cblas.CblasNoTrans;

            const LDA = if (!AH) k else m;
            const LDB = if (!BH) n else k;

            cblas.cblas_zgemm(cblas.CblasRowMajor, ah, bh, m, n, k, &std.math.Complex(T).init(1.0, 0.0), &A.data[0], LDA, &B.data[0], LDB, &std.math.Complex(T).init(0.0, 0.0), &C.data[0], n);
        }
    };
}

/// Eigen struct.
pub fn Eigen(comptime T: type) type {
    return struct {

        /// Logarithm of a matrix.
        pub fn logm(B: *Matrix(T), A: Matrix(T)) void {
            eigen.logm(&B.data[0], &A.data[0], A.rows);
        }
    };
}

/// Expression struct.
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

/// Lapack struct.
pub fn Lapack(comptime T: type) type {
    return struct {

        /// Calculate the condition number of a matrix A, provided its LU decomposition and 1-norm. The result returned.
        pub fn dgecon(ALU: Matrix(T), onorm: T) !T {
            const n: i32 = @intCast(ALU.rows); var rcond: T = 0;

            const info = lapacke.LAPACKE_dgecon(lapacke.LAPACK_ROW_MAJOR, '1', n, &ALU.data[0], n, onorm, &rcond);

            return if (info == 0) rcond else error.ErrorInConditionNumberCalculation;
        }

        /// Calculate the Schur decomposition of matrix A. The JR and JI will contain real and imaginary parts of calculated eigenvalues.
        pub fn dgees(D: *Matrix(T), Q: *Matrix(T), A: Matrix(T), JR: *Matrix(T), JI: *Matrix(T)) !void {
            const n: i32 = @intCast(A.rows); var sdim: i32 = undefined; A.memcpy(D.*);

            const info = lapacke.LAPACKE_dgees(lapacke.LAPACK_ROW_MAJOR, 'V', 'N', null, n, &D.data[0], n, &sdim, &JR.data[0], &JI.data[0], &Q.data[0], n);

            for (0..JR.rows) |i| std.mem.swap(T, JR.ptr(i, i), JR.ptr(0, i));
            for (0..JI.rows) |i| std.mem.swap(T, JI.ptr(i, i), JI.ptr(0, i));

            if (info != 0) return error.ErrorInRealSchurDecomposition;
        }

        /// Solve the linear system Ax = b using the LU decomposition of A. The result is stored in x. The LU decomposition is stored in ALU and the pivot indices are stored in p.
        pub fn dgesv(x: *Vector(T), ALU: *Matrix(T), p: *Vector(i32), A: Matrix(T), b: Vector(T)) !void {
            const n: i32 = @intCast(A.rows); A.memcpy(ALU.*); b.memcpy(x.*);

            const info = lapacke.LAPACKE_dgesv(lapacke.LAPACK_ROW_MAJOR, n, 1, &ALU.data[0], n, &p.data[0], &x.data[0], 1);

            if (info != 0) return error.ErrorInLinearSystemSolution;
        }

        /// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
        pub fn dsyevd(J: *Matrix(T), C: *Matrix(T), A: Matrix(T)) !void {
            const n: i32 = @intCast(A.rows); A.memcpy(C.*); J.fill(0);

            const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

            for (0..J.rows) |i| std.mem.swap(T, J.ptr(i, i), J.ptr(0, i));

            if (info != 0) return error.ErrorInEigenvalueCalculation;
        }

        /// Finds the eigenvalues and eigenvectors of a symmetric-definite generalized eigenproblem A*x = λ*B*x. The eigenvalues are stored in J and the eigenvectors are stored in C. The upper triangular part of Cholesky decmposition will be stored in CH.
        pub fn dsygvd(J: *Matrix(T), C: *Matrix(T), A: Matrix(T), B: Matrix(T), CH: *Matrix(T)) !void {
            const n: i32 = @intCast(A.rows); A.memcpy(C.*); B.memcpy(CH.*); J.fill(0);

            const info = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, &C.data[0], n, &CH.data[0], n, &J.data[0]);

            for (0..J.rows) |i| std.mem.swap(T, J.ptr(i, i), J.ptr(0, i));

            if (info != 0) return error.ErrorInEigenvalueCalculation;
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
